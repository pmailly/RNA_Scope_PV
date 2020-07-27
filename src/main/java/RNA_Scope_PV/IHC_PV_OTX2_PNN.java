package RNA_Scope_PV;


import static Tools.RNAScope_Tools3D.closeImages;
import static Tools.RNAScope_Tools3D.createDonutPop;
import static Tools.RNAScope_Tools3D.dialog;
import static Tools.RNAScope_Tools3D.filterCells;
import static Tools.RNAScope_Tools3D.findAssociatedCell;
import static Tools.RNAScope_Tools3D.findCells;
import static Tools.RNAScope_Tools3D.findPNNCells;
import static Tools.RNAScope_Tools3D.find_background;
import static Tools.RNAScope_Tools3D.readXML;
import static Tools.RNAScope_Tools3D.saveIHCObjects;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
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
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.util.ArrayList;
import javax.xml.parsers.ParserConfigurationException;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;


/*
 * Find nucleus with reference gene and/or virus 
 * based on objects segmentation of gene and virus channels
 * find mRNA in nucleus+ population
 * 
 */

/**
 *
 * @author phm
 */
public class IHC_PV_OTX2_PNN implements PlugIn {
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  static String outDirResults = "";
    public  static String rootName = "";
    public static Calibration cal = new Calibration();
    // threshold to keep PV and Otx2 cells
    public static double PVMinInt, Otx2MinInt;
    public static double sphCell = 0.5;
    public static double minCellVolOtx2 = 600;
    public static double maxCellVolOtx2 = 4000;
    public static double minCellVolPV = 600;
    public static double maxCellVolPV = 8000;
    public static BufferedWriter PV_Analyze, Otx2_Analyze, PNN_Analyze;

    
    
    /** initialize result files
     * 
     */
    private void writeHeaders() throws IOException {        

                // IHC PV results
                FileWriter fwPV = new FileWriter(outDirResults + "PV_results.xls",false);
                PV_Analyze = new BufferedWriter(fwPV);
                // write results headers
                PV_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tPV Mean Intensity\tPV Integrated intensity\tPV Mean background Int\t"
                        + "Std backgroun Int\tPV Corrected Integrated intensity\tOtx2 Corrected Integrated intensity\tPNN Corrected Integrated intensity\n");
                PV_Analyze.flush();

                // IHC Otx2 results
                FileWriter fwOtx2 = new FileWriter(outDirResults + "Otx2_results.xls",false);
                Otx2_Analyze = new BufferedWriter(fwOtx2);
                // write results headers
                Otx2_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tOtx2 Integrated intensity\tOtx2 Mean background Int\t"
                        + "Std background Int\tOtx2 Corrected Integrated intensity\tPV Corrected Integrated intensity\tPNN Corrected Integrated intensity\n");
                Otx2_Analyze.flush();

                 // IHC PNN results
                FileWriter fwPNN = new FileWriter(outDirResults + "PNN_results.xls",false);
                PNN_Analyze = new BufferedWriter(fwPNN);
                // write results headers
                PNN_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tPNN Integrated intensity\tPNN Mean background Int\t"
                        + "Std background Int\tPNN Corrected Integrated intensity\t#PV Cell\tPV Corrected Integrated intensity\t#Otx2 Cell\tOtx2 Corrected Integrated intensity\n");
                PNN_Analyze.flush();
    }
    
    
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing lif Files...");
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
                // Find nd files
                String fileExt = FilenameUtils.getExtension(f);
                
                if (fileExt.equals("nd")) {
                    imageNum++;
                    rootName = FilenameUtils.getBaseName(f);
                    String imageName = inDir + File.separator+f;
                    reader.setId(imageName);
                    // find channel names
                    int chs = reader.getSizeC();
                    String[] ch = new String[chs];
                    if (chs > 1) {
                        for (int n = 0; n < chs; n++)
                            ch[n] = meta.getChannelName(0, n);
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
                 
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(imageName);
                    options.setSplitChannels(true);

                        
                    /** 
                     * read nd
                     * Detect IHC PV cells, measure intensity in PV channel2 and Otx2 channel1
                     * compute donut PV Object and measure in PNN channel0
                     * Detect Otx2 cells measure intensity in Otx2 and PV channels
                     * Detect PNN cells measure intensity in PNN channel and find corresponding PV Cell
                    */

                       // Find xml points file
                        String xmlFile = inDir+ File.separator + rootName + ".xml";
                        String roiFile = inDir+ File.separator + rootName + ".zip";
                        if (!new File(xmlFile).exists() || !new File(roiFile).exists()) {
                            IJ.showStatus("No XML or roi file found !") ;
                        }
                        else {
                            options.setCBegin(0, 0);
                            options.setCEnd(0, 2);
                            
                            // Roi
                            RoiManager rm = new RoiManager(false);
                            rm.runCommand("Open", roiFile);
                            Roi[] rois = rm.getRoisAsArray();
                            // PNN
                            System.out.println("Opening PNN channel ...");
                            ImagePlus imgPNNOrg = BF.openImagePlus(options)[1];
                            //PV
                            System.out.println("Opening PV channel ...");
                            ImagePlus imgPVOrg = BF.openImagePlus(options)[2];
                            //Otx2
                            System.out.println("Opening Otx2 channel ...");
                            ImagePlus imgOtx2Org = BF.openImagePlus(options)[0];
                            // for all rois
                            for (Roi roi : rois) {
                                imgPNNOrg.setRoi(roi);
                                String roiName = roi.getName();
                                ImagePlus imgPNN = new Duplicator().run(imgPNNOrg);
                                // PNN background
                                System.out.println("PNN");
                                double[] bgPNN = find_background(imgPNN);
                                // Find PNN cells with xml points file
                                ArrayList<Point3D> PNNPoints = readXML(xmlFile);
                                Objects3DPopulation PNNPop = findPNNCells(imgPNN, roi, PNNPoints);
                                System.out.println("PNN Cells found : " + PNNPop.getNbObjects());


                                // PV
                                imgPVOrg.setRoi(roi);
                                ImagePlus imgPV = new Duplicator().run(imgPVOrg);
                                //section volume in mm^3
                                double sectionVol = (imgPV.getWidth() * cal.pixelWidth * imgPV.getHeight() * cal.pixelHeight * imgPV.getNSlices() * cal.pixelDepth)/1e9;
                                // PV background
                                System.out.println("PV");
                                double[] bgPV = find_background(imgPV);
                                // find PV cells                          
                                Objects3DPopulation PVPop = findCells(imgPV, roi, 18, 20, 1, "MeanPlusStdDev", true, minCellVolPV, maxCellVolPV);
                                System.out.println("PV Cells found : " + PVPop.getNbObjects());

                                // Otx2
                                imgOtx2Org.setRoi(roi);
                                ImagePlus imgOtx2 = new Duplicator().run(imgOtx2Org);
                                System.out.println("Otx2");
                                // Otx2 background
                                double[] bgOtx2 = find_background(imgOtx2);
                                // Find Otx2 cells
                                Objects3DPopulation Otx2Pop = findCells(imgOtx2, roi, 18, 20, 1, "Huang", true, minCellVolOtx2, maxCellVolOtx2);
                                filterCells(Otx2Pop, 0.55);
                                System.out.println("Otx2 Cells found : " + Otx2Pop.getNbObjects());

                                // save image for objects population
                                saveIHCObjects(PVPop, Otx2Pop, PNNPop, imgPV, outDirResults+rootName+"-"+roiName+"_IHCObjects.tif");    

                                // Compute parameters

                                // PV
                                // create donut
                                float dilatedStepXY = (float) (6/cal.pixelWidth);
                                float dilatedStepZ = (float) (6/cal.pixelDepth);
                                Objects3DPopulation PVDonutPop  = createDonutPop(PVPop, imgPV, dilatedStepXY, dilatedStepZ);
                                ImageHandler imhPV = ImageHandler.wrap(imgPV);
                                ImageHandler imhOtx2 = ImageHandler.wrap(imgOtx2);
                                ImageHandler imhPNN = ImageHandler.wrap(imgPNN);
                                for (int o = 0; o < PVPop.getNbObjects(); o++) {
                                    Object3D obj = PVPop.getObject(o);
                                    Object3D objDonut = PVDonutPop.getObject(o);
                                    double objVol = obj.getVolumeUnit();
                                    double objIntPV = obj.getIntegratedDensity(imhPV);
                                    double objMeanPV = obj.getPixMeanValue(imhPV);
                                    double objIntOtx2 = obj.getIntegratedDensity(imhOtx2);
                                    double objIntPNN = objDonut.getIntegratedDensity(imhPNN);
                                    PV_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+PVPop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objMeanPV+"\t"+objIntPV+"\t"+
                                            bgPV[0]+"\t"+ bgPV[1] + "\t" + (objIntPV - (bgPV[0] * objVol))+"\t"+(objIntOtx2 - (bgOtx2[0] * objVol))+"\t"+(objIntPNN - (bgPNN[0] * objDonut.getVolumeUnit()))+"\n");
                                    PV_Analyze.flush();
                                }

                                // Otx2
                                Objects3DPopulation Otx2DonutPop  = createDonutPop(Otx2Pop, imgOtx2, dilatedStepXY, dilatedStepZ);
                                for (int o = 0; o < Otx2Pop.getNbObjects(); o++) {
                                    Object3D obj = Otx2Pop.getObject(o);
                                    Object3D objDonut = Otx2DonutPop.getObject(o);
                                    double objVol = obj.getVolumeUnit();
                                    double objIntPV = obj.getIntegratedDensity(imhPV);
                                    double objIntOtx2 = obj.getIntegratedDensity(imhOtx2);
                                    double objIntPNN = objDonut.getIntegratedDensity(imhPNN);
                                    Otx2_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+Otx2Pop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objIntOtx2+"\t"+
                                            bgOtx2[0]+"\t"+bgOtx2[1]+"\t"+(objIntOtx2 - (bgOtx2[0] * objVol))+"\t"+(objIntPV - (bgPV[0] * objVol))+"\t"+(objIntPNN - (bgPNN[0] * objDonut.getVolumeUnit()))+"\n");
                                    Otx2_Analyze.flush();
                                }

                                // PNN
                                for (int o = 0; o < PNNPop.getNbObjects(); o++) {
                                    Object3D obj = PNNPop.getObject(o);
                                    double objVol = obj.getVolumeUnit();
                                    double objIntPNN = obj.getIntegratedDensity(imhPNN);
                                    // find associated pv cell
                                    Object3D pvCell = findAssociatedCell(PVPop, obj);
                                    // find associated Otx2 cell
                                    Object3D Otx2Cell = findAssociatedCell(Otx2Pop, obj);
                                    double objIntPV = 0;
                                    double objIntOtx2 = 0;
                                    int pvIndex = -1;
                                    int Otx2Index = -1;
                                    if (pvCell != null) {
                                        objIntPV = pvCell.getIntegratedDensity(imhPV) - (bgPV[0] * pvCell.getVolumeUnit());
                                        pvIndex = PVPop.getIndexOf(pvCell);
                                    }    
                                    if (Otx2Cell != null) {
                                        objIntOtx2 = Otx2Cell.getIntegratedDensity(imhOtx2) - bgOtx2[0] * (Otx2Cell.getVolumeUnit());
                                        Otx2Index = Otx2Pop.getIndexOf(Otx2Cell);
                                    }
                                    PNN_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+PNNPop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objIntPNN+"\t"+
                                            bgPNN[0]+"\t"+bgPNN[1]+"\t"+(objIntPNN - bgPNN[0] * objVol)+"\t"+pvIndex+"\t"+objIntPV+
                                            "\t"+Otx2Index+"\t"+objIntOtx2+"\n");
                                    PNN_Analyze.flush();
                                }
                                closeImages(imgPNN);
                                closeImages(imgOtx2);
                                closeImages(imgPV);
                                
                            }
                            closeImages(imgPVOrg);
                            closeImages(imgOtx2Org);
                            closeImages(imgPNNOrg);
                        }
                    }
                }
                PV_Analyze.close();
                Otx2_Analyze.close();
                PNN_Analyze.close();
            
            } catch (IOException | DependencyException | ServiceException | FormatException | ParserConfigurationException | SAXException ex) {
                Logger.getLogger(IHC_PV_OTX2_PNN.class.getName()).log(Level.SEVERE, null, ex);
            }
        IJ.showStatus("Process done ...");
    }
}
