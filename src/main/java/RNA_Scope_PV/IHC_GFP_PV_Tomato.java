package RNA_Scope_PV;


import static Tools.RNAScope_Tools3D.closeImages;
import static Tools.RNAScope_Tools3D.dialog;
import static Tools.RNAScope_Tools3D.filterCells;
import static Tools.RNAScope_Tools3D.findAssociatedCell;
import static Tools.RNAScope_Tools3D.findCells;
import static Tools.RNAScope_Tools3D.findChannels;
import static Tools.RNAScope_Tools3D.findImageCalib;
import static Tools.RNAScope_Tools3D.findImages;
import static Tools.RNAScope_Tools3D.find_background;
import static Tools.RNAScope_Tools3D.maxCellVol;
import static Tools.RNAScope_Tools3D.minCellVol;
import static Tools.RNAScope_Tools3D.saveIHCObjects;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
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
import java.util.List;
import loci.common.Region;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/*
 * Find PV, OTX2 and PNN cells
 * 
 */

/**
 *
 * @author phm
 */
public class IHC_GFP_PV_Tomato implements PlugIn {
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  static String outDirResults = "";
    public  static String rootName = "";
    public static double sphCell = 0.5;
    public static BufferedWriter Tomato_Analyze;

    
    
    /** initialize result files
     * 
     */
    private void writeHeaders() throws IOException {        

            // IHC Tomato results
            FileWriter fwTomato = new FileWriter(outDirResults + "Tomato_results.xls",false);
            Tomato_Analyze = new BufferedWriter(fwTomato);
            // write results headers
            Tomato_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tTomato cell integrated intensity\tTomato cell Mean background Int\t"
                    + "Tomato cell Std background Int\tTomato cell corrected integrated intensity\tTomato Cell corrected integrated intensity in PV channel\tTomato Cell corrected integrated intensity in GFP channel\t"
                    + "PV cell index\tPV cell corrected integrated intensity in PV channel\tGFP cell index\tGFP cell corrected integrated intensity in GFP channel\n");
            Tomato_Analyze.flush();

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
            
            // Find images with nd extension
            ArrayList<String> imageFile = findImages(imageDir, "nd");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
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
            reader.setId(imageFile.get(0));
            
            // Find channel names
            List<String> channels = findChannels(imageFile.get(0));
            
            // Find image calibration
            Calibration cal = findImageCalib(meta);

            // write headers
            writeHeaders();
            
            // Channels dialog
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("GFP");
            channelsName.add("PV");
            channelsName.add("Tomato");
            if (channels.size() > 1) {
                chs = dialog(channels, channelsName);
                if ( chs == null) {
                    IJ.showStatus("Plugin cancelled");
                    return;
                }
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                ImporterOptions options = new ImporterOptions();
                    
                        
                /** 
                 * read nd
                 * Detect IHC Tomato cells channel3, measure intensity in GFP channel1 and PV channel2
                */

                   // Find roi file
                    String rootFilename =  imageDir+ File.separator + rootName;
                    String roiFile = new File(rootFilename + ".zip").exists() ? rootFilename + ".zip" : rootFilename + ".roi";

                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setQuiet(true);
                    options.setCrop(true);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);


                    // Roi
                    RoiManager rm = new RoiManager(false);
                    if (new File(roiFile).exists())
                        rm.runCommand("Open", roiFile);
                    else
                        rm.add(new Roi(0, 0, reader.getSizeX(), reader.getSizeY()),0);

                    // for all rois
                    for (int r = 0; r < rm.getCount(); r++) {
                        Roi roi = rm.getRoi(r);
                        String roiName = roi.getName();
                        // Find crop region
                        Rectangle rect = roi.getBounds();
                        Region reg = new Region(rect.x, rect.y, rect.width, rect.height);
                        options.setCropRegion(0, reg);

                        roi.setLocation(0, 0);
                        System.out.println("ROI : "+roiName);

                        // Tomato cells
                        System.out.println("Opening Tomato channel " + channels.get(2) +" ...");
                        int channel = channels.indexOf(chs.get(2));
                        ImagePlus imgTomato = BF.openImagePlus(options)[channel];
                        // Tomato background
                        double[] bgTomato = find_background(imgTomato);
                        //section volume in mm^3
                        double sectionVol = (imgTomato.getWidth() * cal.pixelWidth * imgTomato.getHeight() * cal.pixelHeight * imgTomato.getNSlices() * cal.pixelDepth)/1e9;
                        // Find Tomato cells
                        Objects3DPopulation TomatoPop = findCells(imgTomato, roi, 18, 20, 1, "Otsu", false, 0, minCellVol, maxCellVol);
                        filterCells(TomatoPop, 0.55);
                        System.out.println("Tomato Cells found : " + TomatoPop.getNbObjects()  + " in " + roiName);

                        // PV Cells
                        System.out.println("Opening PV channel " + channels.get(1)+ " ...");
                        channel = channels.indexOf(chs.get(1));
                        ImagePlus imgPV = BF.openImagePlus(options)[channel];
                        // PV background
                        double[] bgPV = find_background(imgPV);
                        // find PV cells                          
                        Objects3DPopulation PVPop = findCells(imgPV, roi, 18, 20, 1, "MeanPlusStdDev", false, 0, minCellVol, maxCellVol);
                        System.out.println("PV Cells found : " + PVPop.getNbObjects() + " in " + roiName);

                        // GFP cells
                        System.out.println("Opening GFP channel " + channels.get(0) +" ...");
                        channel = channels.indexOf(chs.get(0));
                        ImagePlus imgGFP = BF.openImagePlus(options)[channel];
                        // GFP background
                        double[] bgGFP = find_background(imgGFP);
                        Objects3DPopulation GFPPop = findCells(imgGFP, roi, 18, 20, 1, "MeanPlusStdDev", true, 20, minCellVol, maxCellVol);
                        System.out.println("GFP Cells found : " + GFPPop.getNbObjects() + " in " + roiName);



                        // save image for objects population
                        saveIHCObjects(GFPPop, PVPop, TomatoPop, imgPV, outDirResults+rootName+"-"+roiName+"_IHCObjects.tif");    

                        // Compute parameters

                        ImageHandler imhPV = ImageHandler.wrap(imgPV);
                        ImageHandler imhTomato = ImageHandler.wrap(imgTomato);
                        ImageHandler imhGFP = ImageHandler.wrap(imgGFP);

                        // Tomato
                        for (int o = 0; o < TomatoPop.getNbObjects(); o++) {
                            Object3D tomatoCell = TomatoPop.getObject(o);
                            double tomatoCellVol = tomatoCell.getVolumeUnit();
                            // find associated PV cell and integrated intensity
                            Object3D pvCell = findAssociatedCell(PVPop, tomatoCell);
                            int pvCellIndex = -1;
                            double pvCellIntChPVCor = 0;
                            if (pvCell != null) {
                                pvCellIntChPVCor = pvCell.getIntegratedDensity(imhPV) - (bgPV[0] * pvCell.getVolumePixels());
                                pvCellIndex = PVPop.getIndexOf(pvCell);
                            }
                            // find associated GFP cell and integrated intensity
                            Object3D GFPCell = findAssociatedCell(GFPPop, tomatoCell);
                            int GFPCellIndex = -1;
                            double GFPCellIntChGFPCor = 0;

                            if (GFPCell != null) {
                                GFPCellIntChGFPCor = GFPCell.getIntegratedDensity(imhGFP) - (bgGFP[0] * GFPCell.getVolumePixels());
                                GFPCellIndex = GFPPop.getIndexOf(GFPCell);
                            }
                            // Find tomato integrated intensity in PV and GFP channel
                            double tomatoCellIntChTomato = tomatoCell.getIntegratedDensity(imhTomato);
                            double tomatoCellIntChTomatoCor = tomatoCellIntChTomato - (bgTomato[0] * tomatoCell.getVolumePixels());
                            double tomatoCellIntChPVCor = tomatoCell.getIntegratedDensity(imhPV) - (bgPV[0] * tomatoCell.getVolumePixels());
                            double tomatoCellIntChGFPCor = tomatoCell.getIntegratedDensity(imhGFP) - (bgGFP[0] * tomatoCell.getVolumePixels());


                            // Write results
                            Tomato_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+TomatoPop.getNbObjects()/sectionVol+"\t"+o+"\t"+tomatoCellVol+"\t"+tomatoCellIntChTomato+"\t"+
                                    bgTomato[0]+"\t"+bgTomato[1]+"\t"+tomatoCellIntChTomatoCor+"\t"+tomatoCellIntChPVCor+"\t"+tomatoCellIntChGFPCor+"\t"+pvCellIndex+"\t"+pvCellIntChPVCor+"\t"+
                                    GFPCellIndex+"\t"+ GFPCellIntChGFPCor+"\n");
                            Tomato_Analyze.flush();
                        }
                        closeImages(imgGFP);
                        closeImages(imgTomato);
                        closeImages(imgPV);

                    }
                }
                if (Tomato_Analyze != null) 
                    Tomato_Analyze.close();
            
            } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
                Logger.getLogger(IHC_GFP_PV_Tomato.class.getName()).log(Level.SEVERE, null, ex);
            }
        IJ.showStatus("Process done ...");
    }
}
