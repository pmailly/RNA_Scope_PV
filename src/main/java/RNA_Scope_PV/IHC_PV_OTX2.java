package RNA_Scope_PV;


import static Tools.RNAScope_Tools3D.closeImages;
import static Tools.RNAScope_Tools3D.dialog;
import static Tools.RNAScope_Tools3D.filterCells;
import static Tools.RNAScope_Tools3D.findCells;
import static Tools.RNAScope_Tools3D.findChannels;
import static Tools.RNAScope_Tools3D.findImageCalib;
import static Tools.RNAScope_Tools3D.findImages;
import static Tools.RNAScope_Tools3D.find_background;
import static Tools.RNAScope_Tools3D.labelsObject;
import static Tools.RNAScope_Tools3D.maxCellVol;
import static Tools.RNAScope_Tools3D.minCellVol;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
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
import ij.plugin.RGBStackMerge;
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
public class IHC_PV_OTX2 implements PlugIn {
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  static String outDirResults = "";
    public  static String rootName = "";
    // threshold to keep PV and Otx2 cells
    public static double PVMinInt, Otx2MinInt;
    public static double sphCell = 0.5;
    public static BufferedWriter PV_Analyze, Otx2_Analyze;

    
    
    /** initialize result files
     * 
     */
    private void writeHeaders() throws IOException {        

                // IHC PV results
                FileWriter fwPV = new FileWriter(outDirResults + "PV_results.xls",false);
                PV_Analyze = new BufferedWriter(fwPV);
                // write results headers
                PV_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tPV Mean Intensity\tPV Integrated intensity\tPV Mean background Int\t"
                        + "Std backgroun Int\tPV Corrected Integrated intensity\tOtx2 Corrected Integrated intensity\n");
                PV_Analyze.flush();

                // IHC Otx2 results
                FileWriter fwOtx2 = new FileWriter(outDirResults + "Otx2_results.xls",false);
                Otx2_Analyze = new BufferedWriter(fwOtx2);
                // write results headers
                Otx2_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tOtx2 Integrated intensity\tOtx2 Mean background Int\t"
                        + "Std background Int\tOtx2 Corrected Integrated intensity\tPV Corrected Integrated intensity\n");
                Otx2_Analyze.flush();
    }
    
    public static void saveIHCObjects(Objects3DPopulation PVPop, Objects3DPopulation Otx2Pop, ImagePlus imgCells, String name) {
       // PVPop blue Otx2Pop/Tomato red PNNPop green
        ImageHandler pvImgObj = ImageHandler.wrap(imgCells).createSameDimensions();
        ImageHandler otx2ImgObj = pvImgObj.duplicate();
        // draw obj population
        if (PVPop != null) {
            PVPop.draw(pvImgObj, 64);
            labelsObject(PVPop, pvImgObj.getImagePlus());
        }
        if (Otx2Pop != null) {
            Otx2Pop.draw(otx2ImgObj, 64);
            labelsObject(Otx2Pop, otx2ImgObj.getImagePlus());
        }
       
        ImagePlus[] imgColors = {otx2ImgObj.getImagePlus(), pvImgObj.getImagePlus(), null, imgCells};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(imgCells.getCalibration());
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        pvImgObj.closeImagePlus(); 
        otx2ImgObj.closeImagePlus();
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
            ArrayList<String> imageFile = findImages(imageDir, "lif");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with lif extension");
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
            
            // write headers
            writeHeaders();
            
            // Find channel names
            List<String> channels = findChannels(imageFile.get(0));
            
            // Find image calibration
            Calibration cal = findImageCalib(meta);

            // write headers
            writeHeaders();
            
            // Channels dialog
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("OTX2");
            channelsName.add("PV");
            
            if (channels.size() > 1) {
                chs = dialog(channels, channelsName);
                if ( chs == null) {
                    IJ.showStatus("Plugin cancelled");
                    return;
                }
            }
            
             /** 
             * 
             * Detect IHC PV cells, measure intensity in PV channel2 and Otx2 channel1
             * Detect Otx2 cells measure intensity in Otx2 and PV channels
            */
             
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);  
                ImporterOptions options = new ImporterOptions();
                String rootFilename =  imageDir+ File.separator + rootName;
                int series=  reader.getSeriesCount();
                // for all series
                for (int s = 0; s < series; s++) {
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setQuiet(true);
                    options.setCrop(true);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    String seriesName = meta.getImageName(s);
                    options.setSeriesOn(s, true);
                    
                    // Roi
                    RoiManager rm = new RoiManager(false);
                    String roiFile = new File(rootFilename + ".zip").exists() ? rootFilename + ".zip" : rootFilename + ".roi";
                    String roiName = "";
                    if (new File(roiFile).exists())
                        rm.runCommand("Open", roiFile);
                    else {
                        rm.add(new Roi(0, 0, reader.getSizeX(), reader.getSizeY()),0);
                        rm.getRoi(0).setName(seriesName);
                    }
                    // for all rois
                    for (int r = 0; r < rm.getCount(); r++) {
                        Roi roi = rm.getRoi(r);
                        roiName = roi.getName();
                        if (roiName.contains(seriesName)) {
                            // Find crop region
                            Rectangle rect = roi.getBounds();
                            Region reg = new Region(rect.x, rect.y, rect.width, rect.height);
                            options.setCropRegion(0, reg);
                            roi.setLocation(0, 0);
                            
                            //PV
                            System.out.println("Opening PV channel ...");
                            int channel = channels.indexOf(chs.get(1));
                            ImagePlus imgPV = BF.openImagePlus(options)[channel];
                            
                            //section volume in mm^3
                            double sectionVol = (imgPV.getWidth() * cal.pixelWidth * imgPV.getHeight() * cal.pixelHeight * imgPV.getNSlices() * cal.pixelDepth)/1e9;
                            // PV background
                            double[] bgPV = find_background(imgPV);
                            // find PV cells                          
                            Objects3DPopulation PVPop = findCells(imgPV, roi, 10, 12, 1, "MeanPlusStdDev", true, 10, minCellVol, maxCellVol);
                            System.out.println("PV Cells found : " + PVPop.getNbObjects() + " in " + roiName);

                            //Otx2
                            System.out.println("Opening Otx2 channel ...");
                            channel = channels.indexOf(chs.get(0));
                            ImagePlus imgOtx2 = BF.openImagePlus(options)[channel];
                            
                            // Otx2 background
                            double[] bgOtx2 = find_background(imgOtx2);
                            // Find Otx2 cells
                            Objects3DPopulation Otx2Pop = findCells(imgOtx2, roi, 10, 12, 1, "Triangle", true, 10, minCellVol, maxCellVol);
                            filterCells(Otx2Pop, 0.55);
                            System.out.println("Otx2 Cells found : " + Otx2Pop.getNbObjects()  + " in " + roiName);

                            // save image for objects population
                            saveIHCObjects(PVPop, Otx2Pop, imgPV, outDirResults+rootName+"-"+roiName+"_IHCObjects.tif");   
                            
                            // Compute parameters

                            // PV
                            ImageHandler imhPV = ImageHandler.wrap(imgPV);
                            ImageHandler imhOtx2 = ImageHandler.wrap(imgOtx2);
                            for (int o = 0; o < PVPop.getNbObjects(); o++) {
                                Object3D obj = PVPop.getObject(o);
                                double objVol = obj.getVolumeUnit();
                                double objIntPV = obj.getIntegratedDensity(imhPV);
                                double objMeanPV = obj.getPixMeanValue(imhPV);
                                double objIntOtx2 = obj.getIntegratedDensity(imhOtx2);
                                PV_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+PVPop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objMeanPV+"\t"+objIntPV+"\t"+
                                        bgPV[0]+"\t"+ bgPV[1] + "\t" + (objIntPV - (bgPV[0] * obj.getVolumePixels()))+"\t"+(objIntOtx2 - (bgOtx2[0] * objVol))+"\n");
                                PV_Analyze.flush();
                            }

                            // Otx2
                            for (int o = 0; o < Otx2Pop.getNbObjects(); o++) {
                                Object3D obj = Otx2Pop.getObject(o);
                                double objVol = obj.getVolumeUnit();
                                double objIntPV = obj.getIntegratedDensity(imhPV);
                                double objIntOtx2 = obj.getIntegratedDensity(imhOtx2);
                                Otx2_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+Otx2Pop.getNbObjects()/sectionVol+"\t"+o+"\t"+objVol+"\t"+objIntOtx2+"\t"+
                                        bgOtx2[0]+"\t"+bgOtx2[1]+"\t"+(objIntOtx2 - (bgOtx2[0] * obj.getVolumePixels()))+"\t"+(objIntPV - (bgPV[0] * obj.getVolumePixels()))+"\n");
                                Otx2_Analyze.flush();
                            }
                            closeImages(imgOtx2);
                            closeImages(imgPV);
                        }
                    }
                    options.setSeriesOn(s, false);
                }
            }
            if (PV_Analyze != null) {
                PV_Analyze.close();
                Otx2_Analyze.close();
            }
            
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(IHC_PV_OTX2.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
