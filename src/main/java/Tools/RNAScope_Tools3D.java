package Tools;


import static RNA_Scope_PV.mRNA_PV.cal;
import static RNA_Scope_PV.mRNA_PV.maxCellVol;
import static RNA_Scope_PV.mRNA_PV.minCellVol;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Measurements;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Object3D_IJUtils;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;



/**
 *
 * @author phm
 */

public class RNAScope_Tools3D {
   
    /**
     *
     * @param img
     */
    public static void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Ask for parameters
     * @param channels
     * @param imageDir
     * @return 
     */
    
    public static ArrayList dialog(String[] channels) {
        ArrayList ch = new ArrayList();
        GenericDialogPlus gd = new GenericDialogPlus("IHC PV parameters");
        gd.addMessage("Choose channels");
        gd.addChoice("Otx2 cells", channels, channels[0]);
        gd.addChoice("PNN cells", channels, channels[1]);
        gd.addChoice("PV cells", channels, channels[2]);
        gd.showDialog();
        ch.add(0, gd.getNextChoice());
        ch.add(1, gd.getNextChoice());
        ch.add(2, gd.getNextChoice());
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
  /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    public static  Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    
    
    /**
     * Find mean Intensity of low PV cells
     * return % of background
     */
    private static double findLowPVIntensity(Objects3DPopulation popPV, double bgInt, ArrayList<Point3D> pvPoints, ImagePlus img) {
        double intensity = 0;
        int nbCell = 0;
        ImageHandler imh = ImageHandler.wrap(img);
        for (int i = 0; i < popPV.getNbObjects(); i++) {
            Object3D obj = popPV.getObject(i);
            // If point inside cell read mean intensity
            for (Point3D p : pvPoints) {
                if (obj.inside(p)) {
                    nbCell++;
                    intensity += obj.getPixMeanValue(imh);
                }
            }
        }
        double bgPourcent = bgInt / (intensity/nbCell);
        System.out.println("PV low = "+ bgPourcent + "of backgroud");
        return(bgPourcent);    
    }
            
            
            
    /**
     * Filters PV cells on sphericity
     */
    public static void filtersPVcells(Objects3DPopulation popPV, double[] bg, ImagePlus img) {
        ImageHandler ima = ImageHandler.wrap(img);
        double bgPV = (bg[0] + bg[1])*2;
        for (int i = 0; i < popPV.getNbObjects(); i++) {
            Object3D obj = popPV.getObject(i);
            double sph = obj.getSphericity(true);
            double intPV = obj.getPixMeanValue(ima);
            if (sph < 0.5 || intPV <= bgPV){
                popPV.removeObject(i);
                i--;
            }
        }
    }
  
    
    
    /**
     * Cells segmentation
     * @param imgCells
     * @param blur1
     * @param blur2
     * @param th
     * @return 
     */
    public static Objects3DPopulation findCells(ImagePlus imgCells, int blur1, int blur2, int med, String th) {
        ImagePlus img = imgCells.duplicate();
        img.setCalibration(imgCells.getCalibration());
        IJ.run(img, "Median...", "radius="+med+" stack");
        ImageStack stack = new ImageStack(imgCells.getWidth(), imgCells.getHeight());
        for (int i = 1; i <= img.getStackSize(); i++) {
            img.setZ(i);
            img.updateAndDraw();
            IJ.run(img, "Nuclei Outline", "blur="+blur1+" blur2="+blur2+" threshold_method="+th+" outlier_radius=15 outlier_threshold=1 max_nucleus_size=100 "
                    + "min_nucleus_size=80 erosion=5 expansion_inner=5 expansion=5 results_overlay");
            img.setZ(1);
            img.updateAndDraw();
            ImagePlus mask = new ImagePlus("mask", img.createRoiMask().getBufferedImage());
            ImageProcessor ip =  mask.getProcessor();
            ip.invertLut();
            for (int n = 0; n < 3; n++) {
                ip.erode();
                ip.dilate();
            }
            stack.addSlice(ip);
        }
        ImagePlus imgStack = new ImagePlus("Nucleus", stack);
        imgStack.setCalibration(imgCells.getCalibration());
        Objects3DPopulation cellPop = new Objects3DPopulation(getPopFromImage(imgStack).getObjectsWithinVolume​(minCellVol, maxCellVol, true));
        cellPop.removeObjectsTouchingBorders(imgStack, false);
        closeImages(imgStack);
        closeImages(img);
        return(cellPop);
    }
   
    /**
     * 
     * @param img
     * @return 
     */
    private static int findMaxZ_intensity(ImagePlus img) {
        double maxMean = 0;
        int z = 0;
        for (int i = 1; i <= img.getNSlices(); i++) {
            img.setSlice(i);
            ImageStatistics meanInt = img.getStatistics(Measurements.MEAN);
            if (meanInt.mean > maxMean) {
                maxMean = meanInt.mean;
                z = i;
            }
        }
        return(z);
    }
    
    /**
     * Crop  cell arround center position
     * Remove X sections (10 microns) arround cell center
    */
    private static ImagePlus cropCell(ImagePlus img, Point3D pt) {
        int crop = Math.round(5.00f/(float)cal.pixelDepth);
        int stackSize = img.getNSlices();
        int ptZ = pt.getRoundZ();
        
        // find start of stack
        int Zstart = 0;
        int start = 0;
        for (int i = 1; i <= crop; i++) {
            Zstart = ptZ - i;
            if (Zstart == 1) {
                start = i;
                break;
            }
        }
        // find stop of stack
        int Zstop = 0;
        int stop = 0;
        for (int i = 1; i <= crop; i++) {
            Zstop = ptZ + i;
            if (Zstop == stackSize) {
                stop = i;
                break;
            }
        }
        int newPtZ = stop - start + 1 -  Math.min(stop, start);
        System.out.println("Crop="+crop+" Stack size ="+stackSize+" Zstart="+Zstart+" Zstop"+Zstop+" pt="+newPtZ);
        // crop
        Rectangle box = new Rectangle(pt.getRoundX() - 100, pt.getRoundY() - 100, 200, 200);
        img.setSlice(newPtZ);
        img.setRoi(box);
        ImagePlus imgCrop = new Duplicator().run(img,Zstart, Zstop);
        imgCrop.setCalibration(img.getCalibration());
        imgCrop.setTitle("cell");
        IJ.run(imgCrop, "Difference of Gaussians", "sigma1=10 sigma2=2 stack");
        Roi roi = new PointRoi(100, 100);
        imgCrop.setSlice(newPtZ);
        imgCrop.setRoi(roi);
        return(imgCrop);
    }
    
    
    /**
     * PNN Cells segmentation
     * find PNN objects inside box 200X200 centrered to point
     * @param imgCells
     * @return 
     */
    public static Objects3DPopulation findPNNCells(ImagePlus imgCells, ArrayList<Point3D> pts) {        
        
        Objects3DPopulation cellPop = new Objects3DPopulation();
        for (int i = 0; i < pts.size(); i++) {
            Point3D pt = pts.get(i);
            ImagePlus imgTmp = cropCell(imgCells, pt);
            IJ.run(imgTmp, "Cell Outliner", "cell_radius=50 tolerance=1 kernel_width=17 kernel_smoothing=1 polygon_smoothing=1 weighting_gamma=3 iterations=3 dilate=10 all_slices");
            ImagePlus cellOutline = WindowManager.getImage("cell Cell Outline");
            cellOutline.hide();
            closeImages(imgTmp);
            Object3D obj = Object3D_IJUtils.createObject3DVoxels(cellOutline, 1);
            obj.setNewCenter(pt.getRoundX(), pt.getRoundY(), pt.getRoundZ());
            cellPop.addObject(obj);
            closeImages(cellOutline);
        }
        return(cellPop);
    }
    
    public static Object3D findAssociatedCell(Objects3DPopulation pop, Object3D cellObj) {
        Object3D objAsso = null;
        for (int i = 0; i < pop.getNbObjects(); i ++) {
            Object3D obj = pop.getObject(i);
            if (cellObj.hasOneVoxelColoc(obj) || cellObj.includesBox(obj))
                objAsso = obj;
        }
        return(objAsso);
    }
    
    public static void saveRNAObjects(Objects3DPopulation cellsPop,ImagePlus imgCells, String name) {
        // green Objects, gray image
        ImageHandler imgObjs = ImageHandler.wrap(imgCells).createSameDimensions();
        // draw obj population
        cellsPop.draw(imgObjs, 64);
        labelsObject(cellsPop, imgObjs.getImagePlus());
        ImagePlus[] imgColors = {null, imgObjs.getImagePlus(), null, imgCells};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        imgObjs.closeImagePlus();
    }
    
    public static void saveIHCObjects(Objects3DPopulation PVPop, Objects3DPopulation TomatoPop, Objects3DPopulation PNNPop, ImagePlus imgCells, String name) {
       // PVPop blue TomatoPop red PNNPop green
        ImageHandler pvImgObj = ImageHandler.wrap(imgCells).createSameDimensions();
        ImageHandler tomatoImgObj = pvImgObj.duplicate();
        ImageHandler pnnImgObj = pvImgObj.duplicate();
        // draw obj population
        if (PVPop != null) {
            PVPop.draw(pvImgObj, 64);
            labelsObject(PVPop, pvImgObj.getImagePlus());
        }
        if (TomatoPop != null)
            TomatoPop.draw(tomatoImgObj, 64);
        if (PNNPop != null) {
            PNNPop.draw(pnnImgObj, 64);
            labelsObject(PNNPop, pnnImgObj.getImagePlus());
        }
        
        ImagePlus[] imgColors = {tomatoImgObj.getImagePlus(), pnnImgObj.getImagePlus(), pvImgObj.getImagePlus(), imgCells};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
        pvImgObj.closeImagePlus(); 
        tomatoImgObj.closeImagePlus();
        pnnImgObj.closeImagePlus();
    }

    
    /**
     * Do Z projection
     * @param img
     * @param projection parameter
     * @return 
     */
    private static ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * @param img 
    */
    public static double[] find_background(ImagePlus img) {
      double[] bg = new double[2];
      img.deleteRoi();
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background =  " + bg[0] + "+-" + bg[1]);
        closeImages(imgProj);
      return(bg);
    }
    
    /**
     * Create donut object population
     * 
     */
    public static Objects3DPopulation createDonutPop(Objects3DPopulation pop, ImagePlus img, float dilateStepXY, float dilateStepZ) {
        ImagePlus imgCopy = img.duplicate();
        ImageInt imgBin = ImageInt.wrap(imgCopy);
        Objects3DPopulation donutPop = new Objects3DPopulation();
        Object3D obj, objDil;
        for (int i = 0; i < pop.getNbObjects(); i++) {
            imgBin.fill(0);
            obj = pop.getObject(i);
            objDil = obj.getDilatedObject(dilateStepXY, dilateStepXY, dilateStepZ);
            objDil.draw(imgBin, 255);
            obj.draw(imgBin, 0);
            Objects3DPopulation tmpPop = getPopFromImage(imgBin.getImagePlus());
            donutPop.addObject(tmpPop.getObject(0));
        }
        closeImages(imgCopy);
        imgBin.closeImagePlus();
        return(donutPop);
    }
    
    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public static void labelsObject (Objects3DPopulation popObj, ImagePlus img) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 24);
        for (int n = 0; n < popObj.getNbObjects(); n++) {
            Object3D obj = popObj.getObject(n);
            int[] box = obj.getBoundingBox();
            int z = (int)obj.getCenterZ();
            int x = box[0] - 2;
            int y = box[2] - 2;
            img.setSlice(z+1);
            ImageProcessor ip = img.getProcessor();
            ip.setFont(tagFont);
            ip.setColor(255);
            ip.drawString(Integer.toString(n), x, y);
            img.updateAndDraw();
        }
    }
    
    /**
     * 
     * @param xmlFile
     * @return
     * @throws ParserConfigurationException
     * @throws SAXException
     * @throws IOException 
     */
    public static ArrayList<Point3D> readXML(String xmlFile) throws ParserConfigurationException, SAXException, IOException {
        ArrayList<Point3D> ptList = new ArrayList<>();
        double x = 0, y = 0 ,z = 0;
        File fXmlFile = new File(xmlFile);
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	Document doc = dBuilder.parse(fXmlFile);
        doc.getDocumentElement().normalize();
        NodeList nList = doc.getElementsByTagName("Marker");
        for (int n = 0; n < nList.getLength(); n++) {
            Node nNode = nList.item(n);
            if (nNode.getNodeType() == Node.ELEMENT_NODE) {
                Element eElement = (Element) nNode;
                x = Double.parseDouble(eElement.getElementsByTagName("MarkerX").item(0).getTextContent());
                y = Double.parseDouble(eElement.getElementsByTagName("MarkerY").item(0).getTextContent());
                z = Double.parseDouble(eElement.getElementsByTagName("MarkerZ").item(0).getTextContent());
            }
            Point3D pt = new Point3D(x, y, z);
            ptList.add(pt);
        }
        return(ptList);
    }
    
}