package fluid;

import processing.core.PApplet;

public class FluidBox extends PApplet {
	
	private static final long serialVersionUID = -6185508289654966583L;
	int width = 500;
	int height = 500;
	int hcells = 100;
	int vcells = 100;
	int pmx;
	int pmy;
	int fmx;
	int fmy;
	int act = 0;
	boolean go = true;
	boolean paint = true;
	float angle = (float) 5.2;
	FluidCom fluid;
	
	public void setup() {
		size(width,height,P2D);
		
		fluid = new FluidCom(hcells, vcells, width, height, this, 500.0f);
	}
	
	public void draw() {
		
		if(go) {	
			float a = 0.001f;//0.6f/ (float) frameRate;
			fluid.update(a);
		}
		
		fluid.draw();
		
		stroke(255,0,0);
		line((float)mouseX, (float)mouseY, (float)mouseX + 30*cos(angle), (float)mouseY + 30*sin(angle));
		stroke(0);
		
		fill(255);
		text((int)frameRate,20,20);
		
		if (mousePressed) {
			if(mouseButton == LEFT && paint) {
				switch(act) {
				case 0:
					//fluid.addDens(mouseX, mouseY);
					fluid.addLinDens(fmx, fmy, mouseX, mouseY);
					break;
				case 1:
					fluid.addTemp(mouseX, mouseY);
					break;
				case 2:
					fluid.addInk(mouseX, mouseY);
					break;
				case 3:
					fluid.addWall(mouseX, mouseY);
					break;
				case 4:
					fluid.remWall(mouseX, mouseY);
					break;
				case 5:
					fluid.printVels(mouseX, mouseY);
					break;
				case 6:
					fluid.addFuel(mouseX, mouseY);
					break;
				}
				
			}
			if((mouseButton == LEFT && !paint) || mouseButton == RIGHT) fluid.addVel(mouseX, mouseY, angle);
			
			fmx = mouseX;
			fmy = mouseY;
		}
		
		pmx = mouseX;
		pmy = mouseY;
	}
	
	public static void main(String args[]) {
		PApplet.main(new String[] { "--present", "fluids.FluidBox" });
	}
	
	/*public void mousePressed() {
		fmx = mouseX;
		fmy = mouseY;
	}*/
	
	public void mouseMoved() {
		fmx = mouseX;
		fmy = mouseY;
	}
	
	public void keyPressed() {
		println("key event detected with num: " + (int)key);
		switch(key) {
		case 32://space
			if(go) go = false;
			else go = true;
			break;
		case 49:
			if(paint) paint = false;
			else paint = true;
			break;
		case 110://n
			angle += 0.2;
			break;
		case 109://m
			angle -= 0.2;
			break;
		case 114://r
			fluid.reset();
			break;
		case 100://d
			act = 2;
			break;
		case 101://e
			act = 6;
			break;
		case 115://s
			act = 1;
			break;
		case 97://a
			act = 0;
			break;
		case 119://w
			act = 3;
			break;
		case 113://q
			act = 4;
			break;
		case 118://v
			fluid.toggleVelocities();
			break;
		case 112://p
			act = 5;
			break;
		case 102://f
			fluid.update(0.1f);
			fluid.debugDraw();
			//fluid.printAccs("Global");
			break;
		case TAB:
			act = (act + 1) % 5;
			break;
		default:
			break;
		}
	}
}
