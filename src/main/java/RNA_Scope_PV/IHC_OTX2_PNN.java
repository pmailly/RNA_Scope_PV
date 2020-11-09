/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package RNA_Scope_PV;

import static Tools.RNAScope_Tools3D.findCells;
import static Tools.RNAScope_Tools3D.closeImages;
import static Tools.RNAScope_Tools3D.createDonutPop;
import static Tools.RNAScope_Tools3D.findAssociatedCell;
import static Tools.RNAScope_Tools3D.findPNNCells;
import static Tools.RNAScope_Tools3D.find_background;
import static Tools.RNAScope_Tools3D.maxCellVol;
import static Tools.RNAScope_Tools3D.median_filter;
import static Tools.RNAScope_Tools3D.minCellVol;
import static Tools.RNAScope_Tools3D.readXML;
import static Tools.RNAScope_Tools3D.saveIHCObjects;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.awt.Rectangle;
import java.util.ArrayList;
import javax.xml.parsers.ParserConfigurationException;
import loci.common.Region;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;
/*
 * Find OTX2 and PNN cells 
 * 
 * @author phm
 */


public class IHC_OTX2_PNN implements PlugIn {
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  static String outDirResults = "";
    public  static String rootName = "";
    public static Calibration cal = new Calibration();
    // threshold to keep OTX2 cells
    public static double OTX2MinInt;
    public static double sphCell = 0.5;
    public static BufferedWriter OTX2_Analyze, PNN_Analyze;
    
    
    /** initialize result files
     * 
     */
    private void writeHeaders() throws IOException {        

                // IHC OTX2 results
                FileWriter fwOTX2 = new FileWriter(outDirResults + "OTX2_results.xls",false);
                OTX2_Analyze = new BufferedWriter(fwOTX2);
                // write results headers
                OTX2_Analyze.write("Image Name\tSeries name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tOTX2 Mean Intensity\tOTX2 Integrated intensity\tOTX2 Mean background Int\t"
                        + "Std backgroun Int\tOTX2 Corrected Integrated intensity\tPNN Corrected Integrated intensity\n");
                OTX2_Analyze.flush();

                 // IHC PNN results
                FileWriter fwPNN = new FileWriter(outDirResults + "PNN_results.xls",false);
                PNN_Analyze = new BufferedWriter(fwPNN);
                // write results headers
                PNN_Analyze.write("Image Name\tSeries Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tPNN Integrated intensity\tPNN Mean background Int\t"
                        + "Std background Int\tPNN Corrected Integrated intensity\t#OTX2 Cell\tOTX2 Corrected Integrated intensity\n");
                PNN_Analyze.flush();
    }
    
    /**
     * Ask for parameters
     * @param channels
     * @param imageDir
     * @return 
     */
    
    public static ArrayList dialog(String[] channels) {
        ArrayList ch = new ArrayList();
        GenericDialogPlus gd = new GenericDialogPlus("IHC OTX2 parameters");
        gd.addMessage("Choose channels");
        gd.addChoice("PNN cells", channels, channels[2]);
        gd.addChoice("OTX2 cells", channels, channels[0]);
        gd.showDialog();
        ch.add(0, gd.getNextChoice());
        ch.add(1, gd.getNextChoice());
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    

    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing Image Files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0;
            ArrayList<String> channels = new ArrayList();
            for (String f : imageFile) {
                // Find lif files
                String fileExt = FilenameUtils.getExtension(f);
                if (fileExt.equals("lif")) {
                    imageNum++;
                    rootName = FilenameUtils.getBaseName(f);
                    String imageName = inDir + File.separator+f;
                    reader.setId(imageName);
                    // find channel names
                    int chs = reader.getSizeC();
                    String[] ch = new String[chs];
                    if (chs > 1) {
                        for (int n = 0; n < chs; n++) 
                            if (meta.getChannelExcitationWavelength(0, n) == null)
                                ch[n] = Integer.toString(n);
                            else 
                                ch[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                    }
                    if (imageNum == 1) {   
                        // read image calibration
                        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
                        cal.pixelHeight = cal.pixelWidth;
                        if (meta.getPixelsPhysicalSizeZ(0) != null)
                            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
                        else
                            cal.pixelDepth = 1;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
                        // dialog
                        if (chs > 1) {
                            channels = dialog(ch);
                            if ( channels == null)
                                return;
                        }
                        
                        // write headers
                        writeHeaders();
                    }
                    
                    // Find roi file name
                    String roiFile = inDir+ File.separator + rootName + ".zip";
                    if (!new File(roiFile).exists()) {
                        IJ.showStatus("No roi file found !") ;
                        return;
                    }
                    else {
                        ImporterOptions options = new ImporterOptions();
                        int series=  reader.getSeriesCount();
                        // for all series
                        for (int s = 0; s < series; s++) {
                            String seriesName = meta.getImageName(s);
                            // Find xml points file
                            String xmlFile = inDir+ File.separator + rootName + "_" + seriesName + ".xml";
                            if (!new File(xmlFile).exists()) {
                                IJ.showStatus("No XML or roi file found !") ;
                            }
                            else {
                                // Roi
                                RoiManager rm = new RoiManager(false);
                                rm.runCommand("Open", roiFile);
                                
                                /** 
                                 * read lif
                                 * Detect IHC OTX2 cells channel 0, measure intensity
                                 * compute donut OTX2 Object and measure in PNN channel 2
                                 * Detect PNN cells measure intensity in PNN channel and find corresponding OTX2 Cell
                                */
                                options.setId(imageName);
                                options.setSplitChannels(true);
                                options.setQuiet(true);
                                options.setCrop(true);
                                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                                options.setCBegin(s, 1);
                                options.setCEnd(s, 2);
                                options.setSeriesOn(s, true);
                                
                                
                                // for all rois
                                for (int r = 0; r < rm.getCount(); r++) {
                                    Roi roi = rm.getRoi(r);
                                    String roiName = roi.getName();
                                    
                                    // if roi name contains series name
                                    // Find crop region
                                    if (roiName.contains(seriesName)) {
                                        Rectangle rect = roi.getBounds();
                                        Region reg = new Region(rect.x, rect.y, rect.width, rect.height);
                                        options.setCropRegion(s, reg);
                                        
                                        
                                        // Find PNN cells with xml points file
                                        ArrayList<Point3D> PNNPoints = readXML(xmlFile, roi);
                                        roi.setLocation(0, 0);
                                        
                                        // PNN
                                        System.out.println("Series : "+seriesName+" ROI : "+roiName);
                                        System.out.println("Opening PNN channel ...");
                                        ImagePlus imgPNN = BF.openImagePlus(options)[1];
                                        // PNN background
                                        double[] bgPNN = find_background(imgPNN);
                                        median_filter(imgPNN, 1);
                                        Objects3DPopulation PNNPop = findPNNCells(imgPNN, roi, PNNPoints);
                                        System.out.println("PNN Cells found : " + PNNPop.getNbObjects());
                                        
                                        //OTX2
                                        System.out.println("Opening OTX2 channel ...");
                                        ImagePlus imgOTX2 = BF.openImagePlus(options)[0];
                                        //section volume in mm^3
                                        double sectionVol = (imgOTX2.getWidth() * cal.pixelWidth * imgOTX2.getHeight() * cal.pixelHeight * imgOTX2.getNSlices() * cal.pixelDepth)/1e9;
                                        // OTX2 background
                                        double[] bgOTX2 = find_background(imgOTX2);
                                        // find OTX2 cells                          
                                        Objects3DPopulation OTX2Pop = findCells(imgOTX2, roi, 4, 6, 1, "MeanPlusStdDev", false, minCellVol, maxCellVol);
                                        System.out.println("OTX2 Cells found : " + OTX2Pop.getNbObjects());
                               
                                        // save image for objects population
                                        saveIHCObjects(OTX2Pop, null, PNNPop, imgOTX2, outDirResults+rootName+"-"+seriesName+"_"+roiName+"_IHCObjects.tif");    
                                        
                                        // Compute parameters

                                        // OTX2
                                        // create donut
                                        float dilatedStepXY = (float) (6/cal.pixelWidth);
                                        float dilatedStepZ = (float) (6/cal.pixelDepth);
                                        Objects3DPopulation OTX2DonutPop  = createDonutPop(OTX2Pop, imgOTX2, dilatedStepXY, dilatedStepZ);
                                        ImageHandler imhOTX2 = ImageHandler.wrap(imgOTX2);
                                        ImageHandler imhPNN = ImageHandler.wrap(imgPNN);
                                        for (int o = 0; o < OTX2Pop.getNbObjects(); o++) {
                                            Object3D obj = OTX2Pop.getObject(o);
                                            Object3D objDonut = OTX2DonutPop.getObject(o);
                                            double objVol = obj.getVolumeUnit();
                                            double objIntOTX2 = obj.getIntegratedDensity(imhOTX2);
                                            double objMeanOTX2 = obj.getPixMeanValue(imhOTX2);
                                            double objIntPNN = objDonut.getIntegratedDensity(imhPNN);
                                            OTX2_Analyze.write(rootName+"\t"+seriesName+"\t"+roiName+"\t"+sectionVol+"\t"+OTX2Pop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objMeanOTX2+"\t"+objIntOTX2+"\t"+
                                                    bgOTX2[0]+"\t"+ bgOTX2[1] + "\t" + (objIntOTX2 - (bgOTX2[0] * obj.getVolumePixels()))+"\t"+(objIntPNN - (bgPNN[0] * objDonut.getVolumePixels()))+"\n");
                                            OTX2_Analyze.flush();
                                        }
                                        
                                         // PNN
                                        for (int o = 0; o < PNNPop.getNbObjects(); o++) {
                                            Object3D obj = PNNPop.getObject(o);
                                            double objVol = obj.getVolumeUnit();
                                            double objIntPNN = obj.getIntegratedDensity(imhPNN);
                                            // find associated pv cell
                                            Object3D pvCell = findAssociatedCell(OTX2Pop, obj);
                                            double objIntOTX2 = 0;
                                            int pvIndex = -1;
                                            if (pvCell != null) {
                                                objIntOTX2 = pvCell.getIntegratedDensity(imhOTX2) - (bgOTX2[0] * pvCell.getVolumePixels());
                                                pvIndex = OTX2Pop.getIndexOf(pvCell);
                                            }    
                                            PNN_Analyze.write(rootName+"\t"+seriesName+"\t"+roiName+"\t"+sectionVol+"\t"+PNNPop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objIntPNN+"\t"+
                                                    bgPNN[0]+"\t"+bgPNN[1]+"\t"+(objIntPNN - bgPNN[0] * obj.getVolumePixels())+"\t"+pvIndex+"\t"+objIntOTX2+"\n");
                                            PNN_Analyze.flush();
                                        }
                                        closeImages(imgPNN);
                                        closeImages(imgOTX2);
                                    }
                                }
                            }
                            options.setSeriesOn(s, false);
                        }
                    }
                }
            }
            if (PNN_Analyze != null) {
                OTX2_Analyze.close();
                PNN_Analyze.close();
            }
        } catch (DependencyException | ServiceException | FormatException | IOException | ParserConfigurationException | SAXException ex) {
            Logger.getLogger(IHC_OTX2_PNN.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
    
}
