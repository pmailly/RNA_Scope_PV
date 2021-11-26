/*
 * #%L
 * Ridge Detection plugin for ImageJ
 * %%
 * Copyright (C) 2014 - 2015 Thorsten Wagner (ImageJ java plugin), 1996-1998 Carsten Steger (original C code), 1999 R. Balasubramanian (detect lines code to incorporate within GRASP)
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */

package ridge;

import java.awt.Color;
import java.awt.Font;
import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.CompositeImage;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;


public class Lines_ implements ExtendedPlugInFilter {
    public double lineWidth = 3.5;
    public double contrastHigh = 230;	
    public double contrastLow = 87;
    public double sigma = 1.51;
    public double lowerThresh = 3.06;
    public double upperThresh = 7.99;
    public double minLength = 0;
    public double maxLength = 0;
    public boolean isDarkLine = false;
    public boolean doCorrectPosition = false;
    public boolean doEstimateWidth = true;
    public boolean doExtendLine = true;
    public boolean showJunctionPoints = false;
    public boolean displayResults = false;
    public boolean addToRoiManager = false;
    public boolean makeBinary = true;
    OverlapOption overlapOption = OverlapOption.NONE;
    public boolean showIDs = false;
    public boolean verbose = false;
    boolean isPreview = false;
    boolean contrastOrLineWidthChangedOnce = false;
    public boolean doStack = true;
    private Options usedOptions = null;
    private Lines_ instance = null;
	
	/** For each frame an ArrayList with the lines of a single frame **/// 
	ArrayList<Lines> result;  
	
	/** For each frame an ArrayList with the junctions of a single frame **/
	ArrayList<Junctions> resultJunction; 
	
	ImagePlus imp;
	
	public Lines_(){
		instance = this;
	}
	
	/**
	 * Access the results from an external plugin
	 * @return
	 */
	public Lines_ getInstance(){
		return instance;
	}
	
	

	
	@Override
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("final")) {
			sortLists();
			//assignLinesToJunctions();
			displayContours();
			if(displayResults){
				createResultsTable(true);
			}
			if(addToRoiManager){
				addToRoiManager();
			}
			if (makeBinary) {
				makeBinary();
			}
			return DONE;

		}
		//Check if the apache commons lang library is available
		try {
		    Class.forName("org.apache.commons.lang3.mutable.MutableLong");
		}
		catch (ClassNotFoundException exception) {
		    IJ.error("Please install apache-commons-lang 3", "It seems that the apache-commons-lang 3 library is not installed on your system. \n "
		    		+ "Download the jar file under https://commons.apache.org/proper/commons-lang/ and copy it to plugins/jars");
		    return DONE;
		}
		
		// Abort macro if no image is open
		if(imp==null){
			IJ.error("No image open");
			return DONE;
		}
		
		/*
		 * Because the optional parameters are easier to determine when the image has an non-inverted LUT, 
		 * images with inverted LUT are reset.
		 */
		if(imp.isInvertedLut()){
			IJ.showMessage("LUT reset", "The LUT of your image is inverted and will now reset for better parameter selection");
			IJ.run(imp, "Invert LUT", "");
		}
		
		Line.resetCounter();
		this.imp = imp;
		result = new ArrayList<Lines>();
		resultJunction = new ArrayList<Junctions>();
		return DOES_8G + DOES_STACKS + FINAL_PROCESSING + PARALLELIZE_STACKS;
	}
	
	private void sortLists(){
		
		Collections.sort(result, new Comparator<Lines>() {

			@Override
			public int compare(Lines o1, Lines o2) {
				// TODO Auto-generated method stub
				if(o1.getFrame()<o2.getFrame())
					return -1;
				if(o1.getFrame()>o2.getFrame())
					return 1;
				return 0;
			}
		});
		Collections.sort(resultJunction,new Comparator<Junctions>() {

			@Override
			public int compare(Junctions o1, Junctions o2) {
				if(o1.getFrame()<o2.getFrame())
					return -1;
				if(o1.getFrame()>o2.getFrame())
					return 1;
				return 0;
			}
		});
	}
	
	public void addToRoiManager(){
		RoiManager rm = RoiManager.getInstance();
		if(rm==null){
			rm = new RoiManager();
			
		}
		for (Lines contours : result) {
			for (Line c : contours) {
			
				float[] x = c.getXCoordinates();
				for(int j = 0; j < x.length; j++){
					x[j] = (float) (x[j] + 0.5);
				}
				float[] y = c.getYCoordinates();
				for(int j = 0; j < y.length; j++){
					y[j] = (float) (y[j] + 0.5);
				}
				
			
				FloatPolygon p = new FloatPolygon(x, y,	c.getNumber());
				Roi r = new PolygonRoi(p, Roi.FREELINE);
				r.setPosition(c.getFrame());
				r.setName("C"+c.getID());
				
				rm.addRoi(r);
				
			}
		}
		for (Junctions junctions : resultJunction) {
			for (Junction j : junctions) {
				
				PointRoi pr = new PointRoi(j.x+0.5,j.y+0.5);
				pr.setName("JP-C"+j.getLine1().getID()+"-C"+j.getLine2().getID());
				pr.setPosition(j.getLine1().getFrame());
				rm.addRoi(pr);
			}
		}
		
		rm.setVisible(true);
		rm.runCommand("UseNames", "true");
		
	}

	private void createResultsTable(boolean showJunctions) {
		ResultsTable rt = ResultsTable.getResultsTable();
		ResultsTable rtSum = new ResultsTable();
		rt.setPrecision(3);
		
		Calibration cal= imp.getCalibration();
		for (Lines contours : result) {
			for (Line c : contours) {
				double meanWidth =0;
				for (int i = 0; i < c.num; i++) {
					rt.incrementCounter();
					rt.addValue("Frame", contours.getFrame());
					
					rt.addValue("Contour ID", c.getID());
					rt.addValue("Pos.", i + 1);
					rt.addValue("X", c.col[i]*cal.pixelWidth);
					
					rt.addValue("Y", c.row[i]*cal.pixelHeight);
					rt.addValue("Length", c.estimateLength()*cal.pixelHeight);
					if(doCorrectPosition && doEstimateWidth){
						rt.addValue("Contrast", Math.abs(c.intensity[i]));
						rt.addValue("Asymmetry", Math.abs(c.asymmetry[i]));
					}
					if(doEstimateWidth){
						
						rt.addValue("Line width", (c.width_l[i]+c.width_r[i])*cal.pixelWidth);
						meanWidth+=c.width_l[i]+c.width_r[i];
						rt.addValue("Angle of normal", c.angle[i]);
					}
					rt.addValue("Class", c.getContourClass().toString().substring(5));
				}
				rtSum.incrementCounter();
				rtSum.addValue("Frame", contours.getFrame());
				rtSum.addValue("Contour ID", c.getID());
				rtSum.addValue("Length", c.estimateLength()*cal.pixelWidth);
				
				if(doEstimateWidth){
					rtSum.addValue("Mean line width",meanWidth/c.num *cal.pixelWidth);
				}
			}
		}

		rt.show("Results");
		rtSum.show("Summary");
		
		if (showJunctions) {
			ResultsTable rt2 = new ResultsTable();
			rt2.setPrecision(0);
			for (Junctions junctions : resultJunction) {
				for (Junction j : junctions) {
					rt2.incrementCounter();
					rt2.addValue("Frame", junctions.getFrame());
					rt2.addValue("Contour ID 1", j.getLine1().getID());//c.get( j.cont1)
					rt2.addValue("Contour ID 2", j.getLine2().getID());
					rt2.addValue("X", j.x*cal.pixelWidth);
					rt2.addValue("Y", j.y*cal.pixelHeight);
				}
			}
			rt2.show("Junctions");
		}
	}
	
	public ImagePlus makeBinary() {
		ImagePlus binary = IJ.createHyperStack(imp.getTitle()+" Detected segments",imp.getWidth(), imp.getHeight(),imp.getNChannels(), imp.getStackSize()/imp.getNChannels(),1,8);
		binary.copyScale(imp);

		ImageProcessor binaryProcessor = binary.getProcessor();
		binaryProcessor.invertLut();
		if (imp.getCompositeMode()>0) {
			((CompositeImage)binary).setLuts(imp.getLuts());
		}
		
		ImageStack is = binary.getImageStack();
		ImageProcessor ip = binary.getProcessor();
		
		for (Lines contours : result) {
			for (Line c : contours) {
				
				float[] x = c.getXCoordinates();
				float[] y = c.getYCoordinates();

				int[] x_poly_r = new int[x.length];
				int[] y_poly_r = new int[x.length];
				
				Polygon LineSurface = new Polygon();

				ip = is.getProcessor(c.getFrame());

				ip.setLineWidth(1);
				ip.setColor(255);
			
				for(int j = 0; j < x.length; j++){
					// this draws the identified line
					if (j >0) {
						ip.drawLine((int) Math.round(x[j-1]), (int) Math.round(y[j-1]),(int) Math.round(x[j]), (int) Math.round(y[j]));
					}
			
					// If Estimate Width is ticked, we also draw the line surface in the binary
					if (doEstimateWidth) {

						double nx = Math.sin(c.angle[j]);
						double ny = Math.cos(c.angle[j]);

						//left point coordinates are directly added to the polygon. right coordinates are saved to be added at the end of the coordinates list
						LineSurface.addPoint((int) Math.round(x[j] - c.width_l[j] * nx),(int) Math.round(y[j] - c.width_l[j] * ny));
						
						x_poly_r[j] = (int) Math.round(x[j] + c.width_r[j] * nx);
						y_poly_r[j] = (int) Math.round(y[j] + c.width_r[j] * ny);
					}
				}
				
				if (doEstimateWidth) {
					// loop to add the right coordinates to the end of the polygon, reversed
					for (int j = 0; j < x.length; j++) {
						if (j < x.length) {
							LineSurface.addPoint(x_poly_r[x.length-1-j],y_poly_r[x.length-1-j]);
						}
					}
					//draw surfaces.
					ip.fillPolygon(LineSurface);
				}		
			}
		}
		if (displayResults)
                    binary.show();
		binary.updateAndDraw();
                return(binary);
	}
		
	private void displayContours() {
		imp.setOverlay(null);
		Overlay ovpoly = new Overlay();
		
		double px, py, nx, ny, px_r = 0, py_r = 0, px_l = 0, py_l = 0;
		double last_w_r, last_w_l;

		// Print contour and boundary
		for (int k = 0; k < result.size(); k++) {
			for (int i = 0; i < result.get(k).size(); i++) {
				FloatPolygon polyMitte = new FloatPolygon();

				FloatPolygon polyR = new FloatPolygon();
				FloatPolygon polyL = new FloatPolygon();
				Line cont = result.get(k).get(i);
				int num_points =  cont.num;
				last_w_r = 0;
				last_w_l = 0;
			
				for (int j = 0; j < num_points; j++) {
					
					px = cont.col[j];
					py = cont.row[j];
					nx = Math.sin(cont.angle[j]);
					ny = Math.cos(cont.angle[j]);
					if (doEstimateWidth) {
						px_r = px + cont.width_r[j] * nx;
						py_r = py + cont.width_r[j] * ny;
						px_l = px - cont.width_l[j] * nx;
						py_l = py - cont.width_l[j] * ny;
					}
					
					polyMitte.addPoint((px + 0.5), (py + 0.5));
					if (doEstimateWidth) {
						if (last_w_r > 0 && cont.width_r[j] > 0) {
							polyR.addPoint((px_r + 0.5), (py_r + 0.5));
						}
						if (last_w_l > 0 && cont.width_l[j] > 0) {
							polyL.addPoint((px_l + 0.5), (py_l + 0.5));
						}
					}
					if (doEstimateWidth) {
						last_w_r = cont.width_r[j];
						last_w_l = cont.width_l[j];
					}
				}
				
				
				PolygonRoi polyRoiMitte = new PolygonRoi(polyMitte,
						Roi.POLYLINE);
				
				polyRoiMitte.setStrokeColor(Color.red);
				int position = result.get(k).getFrame();
				if(!doStack || isPreview){
					position = imp.getCurrentSlice();
				}
			
				polyRoiMitte.setPosition(position);
				ovpoly.add(polyRoiMitte);
				
				
				
				
				if (doEstimateWidth) {
					if(polyL.npoints>1){
						PolygonRoi polyRoiRand1 = new PolygonRoi(polyL,
								Roi.POLYLINE);
						polyRoiRand1.setStrokeColor(Color.green);
						position = result.get(k).getFrame();
						if(!doStack || isPreview){
							position = imp.getCurrentSlice();
						}
						polyRoiRand1.setPosition(position);
						ovpoly.add(polyRoiRand1);
	
						PolygonRoi polyRoiRand2 = new PolygonRoi(polyR,
								Roi.POLYLINE);
						polyRoiRand2.setStrokeColor(Color.green);
						polyRoiRand2.setPosition(position);
						ovpoly.add(polyRoiRand2);
					}
				}
				
				//Show IDs
				if(showIDs){/*
					int posx =  polyMitte.xpoints[0];
					int posy =  polyMitte.ypoints[0];
					if(cont.cont_class == contour_class.cont_start_junc){
						posx =  polyMitte.xpoints[polyMitte.npoints-1];
						posy =  polyMitte.ypoints[polyMitte.npoints-1];
					}
					*/
					
					int posx =  (int)polyMitte.xpoints[polyMitte.npoints/2];
					int posy =  (int)polyMitte.ypoints[polyMitte.npoints/2];
					TextRoi tr = new TextRoi(posx , posy, ""+cont.getID());
					tr.setCurrentFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
					tr.setIgnoreClipRect(true);
					tr.setStrokeColor(Color.orange);
					tr.setPosition(resultJunction.get(k).getFrame());
					ovpoly.add(tr);
				}
			}
		}
		if (showJunctionPoints) {
			// Print junctions
			
			for (int k = 0; k < resultJunction.size(); k++) {
				FloatPolygon pointpoly = new FloatPolygon();
				for (int i = 0; i < resultJunction.get(k).size(); i++) {
					
					pointpoly.addPoint(resultJunction.get(k).get(i).x + 0.5,resultJunction.get(k).get(i).y + 0.5);
				}

				PointRoi pointroi = new PointRoi(pointpoly);
				pointroi.setShowLabels(false);
				int position = resultJunction.get(k).getFrame();
				if(!doStack || isPreview){
					position = imp.getCurrentSlice();
				}
				pointroi.setPosition(position);
				ovpoly.add(pointroi);
			}
		}
		if(ovpoly.size()>0){
			imp.setOverlay(ovpoly);
		}
	}
	
	
	@Override
	public void setNPasses(int nPasses) {
		IJ.showProgress(nPasses, imp.getNSlices());
	}

	@Override
	public void run(ImageProcessor ip) {
		
		if (isPreview) {
			Line.resetCounter();
			result = new ArrayList<Lines>();
			resultJunction = new ArrayList<Junctions>();

		}
		
		LineDetector detect = new LineDetector();
		detect.bechatty = verbose;

		result.add(detect.detectLines(ip, sigma, upperThresh, lowerThresh, minLength,maxLength, isDarkLine, doCorrectPosition, doEstimateWidth, doExtendLine, overlapOption));
		usedOptions = detect.getUsedParamters();
		resultJunction.add(detect.getJunctions());

		if (isPreview) {
			displayContours();
			Line.resetCounter();
			result = new ArrayList<Lines>();
			resultJunction = new ArrayList<Junctions>();
		}
	}
	
	
	

	

	
	/**
	 * Return the detected lines
	 * @return ArrayList of lines for each frame.
	 */
	public ArrayList<Lines> getDetectedLines(){
		return result;
	}
	
	/**
	 * Return the detected junctions
	 * @return ArrayList of junctions for each frame.
	 */
	public ArrayList<Junctions> getDetectedJunctions(){
		return resultJunction;
	}
	/**
	 * Returns the parameter which were used for the analysis
	 * @return The used parameters
	 */
	public Options getParameters(){
		return usedOptions;
	}

    @Override
    public int showDialog(ImagePlus arg0, String arg1, PlugInFilterRunner arg2) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
