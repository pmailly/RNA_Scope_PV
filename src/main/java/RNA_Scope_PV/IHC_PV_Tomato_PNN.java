package RNA_Scope_PV;


import static Tools.RNAScope_Tools3D.closeImages;
import static Tools.RNAScope_Tools3D.createDonutPop;
import static Tools.RNAScope_Tools3D.dialog;
import static Tools.RNAScope_Tools3D.filterCells;
import static Tools.RNAScope_Tools3D.findAssociatedCell;
import static Tools.RNAScope_Tools3D.findCells;
import static Tools.RNAScope_Tools3D.findChannels;
import static Tools.RNAScope_Tools3D.findImageCalib;
import static Tools.RNAScope_Tools3D.findImages;
import static Tools.RNAScope_Tools3D.findPNNCells;
import static Tools.RNAScope_Tools3D.find_background;
import static Tools.RNAScope_Tools3D.maxCellVol;
import static Tools.RNAScope_Tools3D.minCellVol;
import static Tools.RNAScope_Tools3D.readXML;
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
 * Find PV, OTX2 and PNN cells
 * 
 */

/**
 *
 * @author phm
 */
public class IHC_PV_Tomato_PNN implements PlugIn {
    
    private final boolean canceled = false;
    private String imageDir = "";
    public  static String outDirResults = "";
    public  static String rootName = "";
    // threshold to keep PV and Tomato cells
    public static double PVMinInt, TomatoMinInt;
    public static double sphCell = 0.5;
    public static BufferedWriter PV_Analyze, Tomato_Analyze, PNN_Analyze;

    
    
    /** initialize result files
     * 
     */
    private void writeHeaders() throws IOException {        

            // IHC Tomato results
            FileWriter fwTomato = new FileWriter(outDirResults + "Tomato_results.xls",false);
            Tomato_Analyze = new BufferedWriter(fwTomato);
            // write results headers
            Tomato_Analyze.write("Image Name\tLayer name\tSection Volume(mm^3)\tCell density (/mm^3)\t#Cell\tCell Vol\tTomato cell integrated intensity\tTomato cell Mean background Int\t"
                    + "Tomato cell Std background Int\tTomato cell corrected integrated intensity\tTomato Cell corrected integrated intensity in PV channel\tTomato Cell corrected integrated intensity in PNN channel\t"
                    + "PV cell index\tPV cell corrected integrated intensity in PV channel\tPNN cell index\tPNN cell corrected integrated intensity in PNN channel\n");
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
            channelsName.add("Tomato");
            channelsName.add("PNN");
            channelsName.add("PV");
            if (channels.size() > 1) {
                chs = dialog(channels, channelsName);
                if ( chs == null) {
                    IJ.showStatus("Plugin cancelled");
                    return;
                }
            }
            
            // write headers
            writeHeaders();
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                ImporterOptions options = new ImporterOptions();

                /** 
                 * read nd
                 * Detect IHC PV cells, measure intensity in PV channel2 and Tomato channel1
                 * compute donut PV Object and measure in PNN channel0
                 * Detect Tomato cells measure intensity in Tomato and PV channels
                 * Detect PNN cells measure intensity in PNN channel and find corresponding PV Cell
                */

                   // Find xml points file
                    String rootFilename =  imageDir+ File.separator + rootName;
                    String xmlFile = rootFilename + ".xml";
                    String roiFile = new File(rootFilename + ".zip").exists() ? rootFilename + ".zip" : rootFilename + ".roi";
                    if (!new File(xmlFile).exists()) {
                        IJ.showStatus("No XML file found !") ;
                    }
                    else {
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

                            // Find PNN cells with xml points file
                            ArrayList<Point3D> PNNPoints = readXML(xmlFile, roi);
                            roi.setLocation(0, 0);

                            // PNN
                            System.out.println("ROI : "+roiName);
                            System.out.println("Opening PNN channel " + channels.get(1) +" ...");
                            int channel = channels.indexOf(chs.get(1));
                            ImagePlus imgPNN = BF.openImagePlus(options)[channel];
                            // PNN background
                            double[] bgPNN = find_background(imgPNN);
                            Objects3DPopulation PNNPop = findPNNCells(imgPNN, roi, PNNPoints);
                            System.out.println("PNN Cells found : " + PNNPop.getNbObjects() + " in " + roiName);

                            //PV
                            System.out.println("Opening PV channel " + channels.get(2)+ " ...");
                            channel = channels.indexOf(chs.get(2));
                            ImagePlus imgPV = BF.openImagePlus(options)[channel];
                            //section volume in mm^3
                            double sectionVol = (imgPV.getWidth() * cal.pixelWidth * imgPV.getHeight() * cal.pixelHeight * imgPV.getNSlices() * cal.pixelDepth)/1e9;
                            // PV background
                            double[] bgPV = find_background(imgPV);
                            // find PV cells                          
                            Objects3DPopulation PVPop = findCells(imgPV, roi, 18, 20, 1, "MeanPlusStdDev", true, 20, minCellVol, maxCellVol);
                            System.out.println("PV Cells found : " + PVPop.getNbObjects() + " in " + roiName);

                            //Tomato
                            System.out.println("Opening Tomato channel " + channels.get(0) +" ...");
                            channel = channels.indexOf(chs.get(0));
                            ImagePlus imgTomato = BF.openImagePlus(options)[channel];
                            // Tomato background
                            double[] bgTomato = find_background(imgTomato);
                            // Find Tomato cells
                            Objects3DPopulation TomatoPop = findCells(imgTomato, roi, 18, 20, 1, "Triangle", true, 20, minCellVol, maxCellVol);
                            filterCells(TomatoPop, 0.55);
                            System.out.println("Tomato Cells found : " + TomatoPop.getNbObjects()  + " in " + roiName);

                            // save image for objects population
                            saveIHCObjects(PVPop, TomatoPop, PNNPop, imgPV, outDirResults+rootName+"-"+roiName+"_IHCObjects.tif");    

                            // Compute parameters

                            // create donut
                            float dilatedStepXY = (float) (3/cal.pixelWidth);
                            float dilatedStepZ = 3;

                            ImageHandler imhPV = ImageHandler.wrap(imgPV);
                            ImageHandler imhTomato = ImageHandler.wrap(imgTomato);
                            ImageHandler imhPNN = ImageHandler.wrap(imgPNN);

                            // Tomato
                            Objects3DPopulation TomatoDonutPop  = createDonutPop(TomatoPop, imgTomato, dilatedStepXY, dilatedStepZ);
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
                                // find associated PNN cell and integrated intensity
                                Object3D pnnCell = findAssociatedCell(PNNPop, tomatoCell);
                                int pnnCellIndex = -1;
                                double pnnCellIntChPNNCor = 0;

                                if (pnnCell != null) {
                                    pnnCellIntChPNNCor = pnnCell.getIntegratedDensity(imhPNN) - (bgPNN[0] * pnnCell.getVolumePixels());
                                    pnnCellIndex = PNNPop.getIndexOf(pnnCell);
                                }
                                // Find tomato integrated intensity in PV and PNN channel
                                Object3D tomatoCellDonut = TomatoDonutPop.getObject(o);
                                double tomatoCellIntChTomato = tomatoCell.getIntegratedDensity(imhTomato);
                                double tomatoCellIntChTomatoCor = tomatoCellIntChTomato - (bgTomato[0] * tomatoCell.getVolumePixels());
                                double tomatoCellIntChPVCor = tomatoCell.getIntegratedDensity(imhPV) - (bgPV[0] * tomatoCell.getVolumePixels());
                                double tomatoCellIntChPNNCor = tomatoCellDonut.getIntegratedDensity(imhPNN) - (bgPNN[0] * tomatoCellDonut.getVolumePixels());


                                // Write results
                                Tomato_Analyze.write(rootName+"\t"+roiName+"\t"+sectionVol+"\t"+TomatoPop.getNbObjects()/sectionVol+"\t"+o+"\t"+tomatoCellVol+"\t"+tomatoCellIntChTomato+"\t"+
                                        bgTomato[0]+"\t"+bgTomato[1]+"\t"+tomatoCellIntChTomatoCor+"\t"+tomatoCellIntChPVCor+"\t"+tomatoCellIntChPNNCor+"\t"+pvCellIndex+"\t"+pvCellIntChPVCor+"\t"+
                                        pnnCellIndex+"\t"+ pnnCellIntChPNNCor+"\n");
                                Tomato_Analyze.flush();
                            }
                            closeImages(imgPNN);
                            closeImages(imgTomato);
                            closeImages(imgPV);
                        }

                    }
                }
                if (Tomato_Analyze != null) 
                    Tomato_Analyze.close();
            
            } catch (IOException | DependencyException | ServiceException | FormatException | ParserConfigurationException | SAXException ex) {
                Logger.getLogger(IHC_PV_Tomato_PNN.class.getName()).log(Level.SEVERE, null, ex);
            }
        IJ.showStatus("Process done ...");
    }
}
