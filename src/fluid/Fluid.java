package fluid;

import processing.core.*;

public class Fluid {
	PApplet parent;
	private int cellw;
	private int cellh;
	private int rowSize;
	private int colSize;
	private int tsize;
	private float[] dens0;
	private float[] dens1;
	private float[] vx0;
	private float[] vx1;
	private float[] vx2;
	private float[] vy0;
	private float[] vy1;
	private float[] vy2;
	private float[] vxs;
	private float[] vys;
	private float[] req;
	private float[] temp0;
	private float[] temp1;
	private float[] ink0;
	private float[] ink1;
	private float[] indReq;
	private float[] noise;
	private float[] nx;
	private float[] ny;
	private int[] sourceCell;
	private boolean[] wall;
	private int[] numParticles;
	private float baseDens;
	private float addVel;
	private boolean pv = true;
	private float nt = 0; 
	
	public float[] getVx0() {
		return vx0;
	}

	public float[] getVy0() {
		return vy0;
	}
	
	public int[] getNumParticles() {
		return numParticles;
	}
	
	public void toggleVelocities() {
		if(pv) pv = false;
		else pv = true;
	}
	
	public Fluid(int cw, int ch, int w, int h, PApplet p, float bd) {
		rowSize = cw + 2;
		colSize = ch + 2;
		parent = p;
		cellw = Math.round(w / cw);
		cellh = Math.round(h / ch);
		
		tsize = rowSize*colSize;
		dens0 = new float[tsize];
		dens1 = new float[tsize];
		vx0 = new float[tsize];
		vx1 = new float[tsize];
		vx2 = new float[tsize];
		vxs = new float[tsize];
		vy0 = new float[tsize];
		vy1 = new float[tsize];
		vy2 = new float[tsize];
		vys = new float[tsize];
		req = new float[tsize];
		temp0 = new float[tsize];
		temp1 = new float[tsize];
		ink0 = new float[tsize];
		ink1 = new float[tsize];
		noise = new float[tsize];
		nx = new float[tsize];
		ny = new float[tsize];
		indReq = new float[4*tsize];
		sourceCell = new int[tsize];
		wall = new boolean[tsize];
		numParticles = new int[tsize];
		baseDens = bd;
		addVel = 100;
		for(int i=0; i<tsize; i++) dens0[i] = bd;
		
		recalcNoiseFields(noise, nx, ny);
	}
	
	public void update(float a) {
		clearWall();
		cdiffuse(dens0, dens1, (float) (0.1*a));
		cdiffuse(temp0, temp1, (float) (1*a));
		cdiffuse(ink0, ink1, (float) (0.1*a));
		
		wipe(vxs); wipe(vys);
		pressure(dens1, vxs, vys, (float) (0.05*a));
		temp(temp1, vxs, vys, (float) (0.1*a));
		fall(ink1, vxs, vys, (float) (0.00001*a));
		addFields(vx0, vy0, vxs, vys);
		
		cdiffuse(vx0, vx2, (float) (0.01*a));
		cdiffuse(vy0, vy2, (float) (0.01*a));
		
		edgeVelocities(vx2, vy2);
		
		//addTurbulence(vx2, vy2, nx, ny, (float) (0.001*a), (float) 1);
		vconf(vx2, vy2, vx0, vy0, vx1, vy1, (float) (0.1*a)); 
		edgeVelocities(vx2, vy2);
		
		aproject(vx1, vy1, vx0, vy0);
		edgeVelocities(vx1, vy1);
		
		friction(vx2, vy2, vx1, vy1, (float) (10*a), (float) (0.1*a), (float) (1*a));
		edgeVelocities(vx1, vy1);
		
		copy(vx1,vx0);
		copy(vy1,vy0);
				
		advect(vx0, vy0, vx1, vx2, (float) (0.05*a));
		advect(vx0, vy0, vy1, vy2, (float) (0.05*a));
		
		edgeVelocities(vx2, vy2);
		
		copy(vx2, vx1);
		copy(vy2, vy1);
		
		badvect(vx1, vy1, vx1, vx2, (float) (0.05*a));
		badvect(vx1, vy1, vy1, vy2, (float) (0.05*a));
		
		edgeVelocities(vx2, vy2);
		
		copy(vx2, vx0);
		copy(vy2, vy0);
		
		copy(dens1, dens0);
		advect(vx0, vy0, dens1, dens0, (float) (0.05*a)); globalDiff(dens0, dens1, (float) (0.01*a));
		buadvect(vx0, vy0, dens1, dens0, (float) (0.05*a)); globalDiff(dens0, dens1, (float) (0.01*a));
		lbound(dens0, baseDens/2);
		
		copy(temp1, temp0);
		advect(vx0, vy0, temp1, temp0, (float) (0.1*a)); copy(temp0, temp1);
		buadvect(vx0, vy0, temp1, temp0, (float) (0.1*a)); swap(temp0, temp1);
		
		copy(ink1, ink0);
		advect(vx0, vy0, ink1, ink0, (float) (0.1*a)); copy(ink0, ink1);
		buadvect(vx0, vy0, ink1, ink0, (float) (0.1*a)); swap(ink0, ink1);
		
	}
	
	public void draw() {
		//leave 1 cell border around the screen
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				int pxu = cell(x+1,y);
				int pxd = cell(x-1,y);
				int pyu = cell(x,y+1);
				int pyd = cell(x,y-1);
				int dens = (int) (dens0[p] + 0.25*(dens0[pxu] + dens0[pxd] + dens0[pyu] + dens0[pyd]));
				int ink = (int) (ink0[p] + 0.25*(ink0[pxu] + ink0[pxd] + ink0[pyu] + ink0[pyd]));
				int temp = (int) (temp0[p] + 0.25*(temp0[pxu] + temp0[pxd] + temp0[pyu] + temp0[pyd]));
				int wallc = 0;
				if(wall[p]) wallc = 255;
				parent.fill(dens + temp, dens + wallc, dens + ink);
				parent.rect((x-1)*cellw, (y-1)*cellh, cellw, cellh);
				if(pv) {
					if(Math.abs(vx0[p]) > 0.01 || Math.abs(vy0[p]) > 0.01) {
						parent.stroke(0,255,0);
						parent.line((float)((x-1)*cellw+0.5*cellw), 
								(float)((y-1)*cellh+0.5*cellh), 
								(float)((x-1)*cellw+0.5*cellw+5*vx0[p]), 
								(float)((y-1)*cellh+0.5*cellh+5*vy0[p]));
						parent.stroke(0);
					}
				}
			}
		}
	}
	
	public void debugDraw() {
		float cw = cellw * (float)(rowSize - 2)/rowSize;
		float ch = cellh * (float)(colSize - 2)/colSize;
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				int dens = (int) (dens0[p]);
				int ink = (int) (ink0[p]);
				int temp = (int) (temp0[p]);
				int wallc = 0;
				if(wall[p]) wallc = 255;
				parent.fill(dens + temp, dens + wallc, dens + ink);
				parent.rect(x*cw, y*ch, cw, ch);
				if(pv) {
					if(Math.abs(vx0[p]) > 0.01 || Math.abs(vy0[p]) > 0.01) {
						parent.stroke(192, 0, 0);
						parent.line((float)(x*cw+0.5*cw), 
								(float)(y*ch+0.5*ch), 
								(float)(x*cw+0.5*cw+5*vx0[p]), 
								(float)(y*ch+0.5*ch+5*vy0[p]));
						parent.stroke(0);
					}
					/*if(Math.abs(nx[p]) > 0.01 || Math.abs(ny[p]) > 0.01) {
						parent.stroke(0, 128, 0);
						parent.line((float)(x*cw+0.5*cw), 
								(float)(y*ch+0.5*ch), 
								(float)(x*cw+0.5*cw+5*nx[p]), 
								(float)(y*ch+0.5*ch+5*ny[p]));
						parent.stroke(0);
					}*/
				}
			}
		}
	}
	
	public void printVels(int mouseX, int mouseY) {
		int posx = (int)((double)mouseX / 600.0 *62);
		int posy = (int)((double)mouseY / 600.0 *62);
		int p = cell(posx, posy);
		PApplet.println(vx0[p] + " " + vy0[p] + " " + vx1[p] + " " + vy1[p] + " " + vx2[p] + " " + vy2[p]);
		PApplet.println(posx + " , " + posy + " , " +dens0[p]);
	}
	
	public void addDens(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		int dens = (int) dens0[p];
		dens += baseDens*4;
		dens0[p] = (float) dens;
	}
	
	public void addDens(int mouseX, int mouseY, float bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float dens = dens0[p];
		dens += 10*bd;
		dens0[p] = dens;
	}
	
	public void addTemp(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = temp0[p];
		temp += baseDens;
		temp0[p] = temp;
	}
	
	public void addTemp(int mouseX, int mouseY, float bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = temp0[p];
		temp += bd;
		temp0[p] = temp;
	}
	
	public void addInk(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float ink = ink0[p];
		ink += baseDens;
		ink0[p] = ink;
	}
	
	public void addVel(int mouseX, int mouseY, float angle) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float nextX = (float) (addVel*Math.cos(angle));
		float nextY = (float) (addVel*Math.sin(angle));
		vx0[p] = nextX;
		vy0[p] = nextY;
	}
	
	public void addWall(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		wall[p] = true;
	}
	
	public void remWall(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		if(wall[p]) {
			wall[p] = false;
			float acc = 0;
			float div = 0;
			int pxu = ncell(posx, posy, posx+1,posy);
			if(pxu != p) {
				acc += dens0[pxu];
				div += 1.0;
			}
			int pxd = ncell(posx, posy, posx-1,posy);
			if(pxd != p) {
				acc += dens0[pxd];
				div += 1.0;
			}
			int pyu = ncell(posx, posy, posx,posy+1);
			if(pyu != p) {
				acc += dens0[pyu];
				div += 1.0;
			}
			int pyd = ncell(posx, posy, posx,posy-1);
			if(pyd != p) {
				acc += dens0[pyd];
				div += 1.0;
			}
			if(div > 0) dens0[p] = acc / div;
			else dens0[p] = 10;
		}
	}
	
	public void copy(float[] ain, float[] aout) {
		for(int i = 0; i < tsize; i++) {
			aout[i] = ain[i];
		}
	}
	
	public void swap(float[] ain, float[] aout) {
		float[] tmp = aout;
		aout = ain;
		ain = tmp;
	}
	
	public void clearWall() {
		for(int i = 0; i < tsize; i++) {
			if(wall[i]) {
				dens0[i] = 0;
				dens1[i] = 0;
				vx0[i] = 0;
				vx1[i] = 0;
				vx2[i] = 0;
				vy0[i] = 0;
				vy1[i] = 0;
				vy2[i] = 0;
			}
		}
	}
	
	public void cutoff(float[] d) {
		for(int i = 0; i < tsize; i++) {
			if(d[i] < 0.001) d[i] = 0;
		}
	}
	
	public void ubound(float[] d, float bnd) {
		for(int i = 0; i < tsize; i++) {
			if(d[i] > bnd) d[i] = bnd;
		}
	}
	public void lbound(float[] d, float bnd) {
		for(int i = 0; i < tsize; i++) {
			if(d[i] < bnd) d[i] = bnd;
		}
	}
	
	public void wipe(float[] d) {
		for(int i = 0; i < tsize; i++) d[i] = 0;
	}
	
	private void addFields(float[] dx, float[] dy, float[] sx, float[] sy) {
		for(int i = 0; i < tsize; i++) {
			dx[i] += sx[i];
			dy[i] += sy[i];
		}
		
	}
	
	public void reset() {
		for(int i = 0; i < tsize; i++) {
			wall[i] = false;
			dens0[i] = baseDens;
			dens1[i] = baseDens;
			temp0[i] = 0;
			temp1[i] = 0;
			ink0[i] = 0;
			ink1[i] = 0;
			vx0[i] = 0;
			vx1[i] = 0;
			vx2[i] = 0;
			vy0[i] = 0;
			vy1[i] = 0;
			vy2[i] = 0;
		}
		recalcNoiseFields(noise, nx, ny);
	}
	
	private void recalcNoiseFields(float[] an, float[] ax, float[] ay) {
		float nscale = (float) 0.03;
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				an[cell(x,y)] = parent.noise((float)(x*nscale), (float)(y*nscale), nt)*255;
			}
		}
		
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int pn = cell(x, y);
				int pxu = cell(x+1, y);
				int pxd = cell(x-1, y);
				int pyu = cell(x, y+1);
				int pyd = cell(x, y-1);
				ax[pn] = (float) (0.5*(an[pyu] - an[pyd]));
				ay[pn] = (float) (-0.5*(an[pxu] - an[pxd]));
			}
		}
		
		nt += 10*nscale;
	}
	
	public void addTurbulence(float[] vx, float[] vy, float[] nx, float[] ny, float c, float d) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int i = cell(x,y);
				float len = (float) Math.sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
				float scalex = (float) - x * (x - rowSize) / (float) (rowSize*rowSize);
				float scaley = (float) - y * (y - colSize) / (float) (colSize*colSize);
				if (len == 0) continue;
				float sc = len*len*len / (len*len*len + d);
				vx[i] = vx[i]+(float) (c*sc*scalex*scalex*nx[i]);
				vy[i] = vy[i]+(float) (c*sc*scaley*scaley*ny[i]);
			}
		}
		
		//recalcNoiseFields(noise, nx, ny);
	}
	
	public void friction(float[] vxa, float[] vya, float[] vxb, float[] vyb, float a, float b, float c) {
		for(int i = 0; i < tsize; i++) {
			float len = (float) Math.sqrt(vxa[i]*vxa[i]+vya[i]*vya[i]);
			if(len == 0) continue;
			float flen = len - len*len*a - len*b - c;
			if(flen < 0) flen = 0;
			float nlen = flen/len;
			vxb[i] = vxa[i] * nlen;
			vyb[i] = vya[i] * nlen;
		}
	}
	
	public void pressure(float[] dens, float[] vx, float[] vy, float a) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				float ax = Math.abs(dens[pxd]) - Math.abs(dens[pxu]);
				float ay = Math.abs(dens[pyd]) - Math.abs(dens[pyu]);
				vx[p] += a*ax;
				vy[p] += a*ay;
			}
		}
	}
	
	public void pressurebias(float[] dens, float[] vx, float[] vy, float a) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				//int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				//int pyd = ncell(x, y, x,y-1);
				/*float ax = Math.abs(dens[pxd]) - Math.abs(dens[pxu]);
				float ay = Math.abs(dens[pyd]) - Math.abs(dens[pyu]);*/
				float ax = dens[p] - dens[pxu];
				float ay = dens[p] - dens[pyu];
				vx[p] += a*ax;
				vy[p] += a*ay;
				vx[pxu] += a*ax;
				vy[pyu] += a*ay;
			}
		}
	}
	
	public void temp(float[] dens, float[] vx, float[] vy, float a) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				float ay = (float) Math.abs(dens[p]);
				vy[p] -= a*ay;
			}
		}
	}
	
	public void fall(float[] dens, float[] vx, float[] vy, float a) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				float ay = (float) Math.abs(dens[p]);
				vy[p] += a*ay;
			}
		}
	}
	
	public void globalDiff(float[] ain, float[] aout, float k) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				aout[p] = (1-4*k)*ain[p] + k*(baseDens);
			}
		}
	}
	
	public void edgeVelocities(float[] vx, float[] vy) {
		for(int x = 0; x < rowSize; x++) {
			int i = cell(x,0);
			if(vy[i] < 0) vy[i] = -vy[i];
			int j = cell(x,colSize-1);
			if(vy[j] > 0) vy[j] = -vy[j];
		}
		for(int y = 0; y < colSize; y++) {
			int i = cell(0,y);
			if(vx[i] < 0) vx[i] = -vx[i];
			int j = cell(rowSize-1,y);
			if(vx[j] > 0) vx[j] = -vx[j];
		}
		
		/***********************
		 * velocities for walls
		 ***********************/
		
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) {
					int pxu = ncell(x, y, x+1,y);
					int pxd = ncell(x, y, x-1,y);
					int pyu = ncell(x, y, x,y+1);
					int pyd = ncell(x, y, x,y-1);
					if(vx[pxu] < 0 && !wall[pxu]) vx[pxu] = -vx[pxu];
					if(vx[pxd] > 0 && !wall[pxd]) vx[pxd] = -vx[pxd];
					if(vy[pyu] < 0 && !wall[pyu]) vy[pyu] = -vy[pyu];
					if(vy[pyd] > 0 && !wall[pyd]) vy[pyd] = -vy[pyd];
				}
			}
		}
	}
	
	public void vconf(float[] vxa, float[] vya, float[] wx, float[] wy, float[] wa, float[] wu, float c) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				wx[p] = (float) (0.5*(vxa[pyu] - vxa[pyd]));
				wy[p] = (float) (-0.5*(vya[pxu] - vya[pxd]));
				wa[p] = (float) Math.sqrt(wx[p]*wx[p] + wy[p]*wy[p]);
			}
		}
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				float dwx = (float) (0.5*(wa[pxu] - wa[pxd]));
				float dwy = (float) (0.5*(wa[pyu] - wa[pyd]));
				float dn = 0;
				float n = (float) Math.sqrt(dwx*dwx+dwy*dwy);
				if(n != 0) {
					dn = (float) (1.0/n);
				}
				vxa[p] += c*(dwy * dn * wx[p]);
				vya[p] -= c*(dwx * dn * wy[p]);
			}
		}
		
	}
	
	public void cdiffuse(float[] ain, float[] aout, float k) {
			for(int x = 0; x < rowSize; x++) {
				for(int y = 0; y < colSize; y++) {
					int i = cell(x,y);
					if(wall[i]) continue;
					int pxu = ncell(x, y, x+1,y);
					int pxd = ncell(x, y, x-1,y);
					int pyu = ncell(x, y, x,y+1);
					int pyd = ncell(x, y, x,y-1);
					aout[i] = (1-4*k)*ain[i] + k*(ain[pxu] + ain[pxd] + ain[pyu] + ain[pyd]);
				}
			}
	}
	
	public void aproject(float[] vx, float[] vy, float[] div, float[] p) {
		for(int x = 1; x < rowSize-1; x++) {
			for(int y = 1; y < colSize-1; y++) {
				int i = cell(x,y);
				int pxu = cell(x+1, y);
				int pxd = cell(x-1, y);
				int pyu = cell(x, y+1);
				int pyd = cell(x, y-1);
				div[i] = (float) (-0.5*(vx[pxu] - vx[pxd] + vy[pyu] - vy[pyd]));
				p[i] = 0;
			}
		}
		setBnd(p, 0);
		setBnd(div, 0);
		
		for(int n=0; n< 20; n++) {
			for(int x = 1; x < rowSize-1; x++) {
				for(int y = 1; y < colSize-1; y++) {
					int i = cell(x,y);
					int pxu = cell(x+1, y);
					int pxd = cell(x-1, y);
					int pyu = cell(x, y+1);
					int pyd = cell(x, y-1);
					p[i] = (div[i] + p[pxu] + p[pxd] + p[pyu] + p[pyd])/4;
				}
			}
			setBnd(p, 0);
		}
		
		for(int x = 1; x < rowSize-1; x++) {
			for(int y = 1; y < colSize-1; y++) {
				int i = cell(x,y);
				int pxu = cell(x+1, y);
				int pxd = cell(x-1, y);
				int pyu = cell(x, y+1);
				int pyd = cell(x, y-1);
				vx[i] += 0.5*(p[pxu]-p[pxd]);
				vy[i] += 0.5*(p[pyu]-p[pyd]);
			}
		}	
	}
	
	
	public void advect(float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				float vx = vxa[p];
				float vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					float xn = x + vx * c;
					float yn = y + vy * c;
					
					float lb = (float) (rowSize - 1.0001);
					float ub = (float) (colSize - 1.0001);
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					if(wall[pn]) {
						//if we advect into a wall, we advect instead into ourselves
						xn = x;
						yn = y;
						if(xn < 0) xn = 0;
						else if (xn > lb) xn = lb;
						if(yn < 0) yn = 0;
						else if (yn > ub) yn = ub;
						
						pn = cell((int)xn,(int)yn);
					}
					
					float fx = xn - (int)xn;
					float fy = yn - (int)yn;
					
					float in = ain[p];
					
					float ia = (float) ((1.0-fy) * (1.0-fx) * in);
					float ib = (float) ((1.0-fy) * fx * in);
					float ic = (float) (fy * (1.0-fx) * in);
					float id = fy * fx * in;
					aout[p] -= (ia+ib+ic+id);
					aout[pn] += ia;
					aout[pn+1] += ib;
					aout[pn+rowSize] += ic;
					aout[pn+rowSize+1] += id;
				}
			}
		}
	}
	
	public void badvect(float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				float vx = vxa[p];
				float vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					float xn = x - vx * c;
					float yn = y - vy * c;
					
					float lb = (float) (rowSize - 1.0001);
					float ub = (float) (colSize - 1.0001);
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					float fx = xn - (int)xn;
					float fy = yn - (int)yn;
					
					float ia = (float) ((1.0-fy) * (1.0-fx) * ain[pn]);
					float ib = (float) ((1.0-fy) * fx * ain[pn+1]);
					float ic = (float) (fy * (1.0-fx) * ain[pn+rowSize]);
					float id = fy * fx * ain[pn+rowSize+1];
					
					aout[p] += (ia+ib+ic+id);
					aout[pn] -= ia;
					aout[pn+1] -= ib;
					aout[pn+rowSize] -= ic;
					aout[pn+rowSize+1] -= id;
				}
			}
		}
	}
	
	public void buadvect(float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		for(int p = 0; p < tsize; p++) req[p]= 0;
		
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				float vx = vxa[p];
				float vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					float xn = x - vx * c;
					float yn = y - vy * c;
					
					float lb = (float) (rowSize - 1.0001);
					float ub = (float) (colSize - 1.0001);
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					float fx = xn - (int)xn;
					float fy = yn - (int)yn;
					
					float ia = (float) ((1.0-fy) * (1.0-fx));
					float ib = (float) ((1.0-fy) * fx);
					float ic = (float) (fy * (1.0-fx));
					float id = fy * fx;
					
					//save sources and what each loses
					sourceCell[p] = pn;
					
					indReq[4*p] = ia;
					indReq[4*p+1] = ib;
					indReq[4*p+2] = ic;
					indReq[4*p+3] = id;
					
					//accumulate how much each cell loses in total
					req[pn] += ia;
					req[pn+1] += ib;
					req[pn+rowSize] += ic;
					req[pn+rowSize+1] += id;
				} else {
					sourceCell[p] = -1;
				}
			}
		}
		
		for(int p = 0; p < tsize; p++) {
			if(wall[p]) continue;
			int pn = sourceCell[p];
			if(pn != -1) {
				//recover previous data
				float ia = indReq[4*p];
				float ib = indReq[4*p+1];
				float ic = indReq[4*p+2];
				float id = indReq[4*p+3];
				
				//get total fractions and rescale requests
				float fa = req[pn];
				float fb = req[pn+1];
				float fc = req[pn+rowSize];
				float fd = req[pn+rowSize+1];
				
				if (fa<1.0f) fa = 1.0f;
				if (fb<1.0f) fb = 1.0f;
				if (fc<1.0f) fc = 1.0f;
				if (fd<1.0f) fd = 1.0f;
				
				ia = ia * ain[pn] / fa;
				ib = ib * ain[pn+1] / fb;
				ic = ic * ain[pn+rowSize] / fc;
				id = id * ain[pn+rowSize+1] / fd;
				
				aout[p] += (ia+ib+ic+id);
				aout[pn] -= ia;
				aout[pn+1] -= ib;
				aout[pn+rowSize] -= ic;
				aout[pn+rowSize+1] -= id;				
			}
		}
		
		for(int p = 0; p < tsize; p++) {
			if(Math.abs(aout[p]) < 1e-8) aout[p] = 0;
		}
	}
	
	private void setBnd(float[] de, int mode) {
		for(int x = 1; x < rowSize - 1; x++) {
			int i = cell(x,0);
			int ip = cell(x,1);
			if(mode == 2) de[i] = -de[ip];
			else de[i] = de[ip];
			int j = cell(x,colSize-1);
			int jp = cell(x,colSize-2);
			if(mode == 2) de[j] = -de[jp];
			else de[j] = de[jp];
		}
		for(int y = 1; y < colSize - 1; y++) {
			int i = cell(0,y);
			int ip = cell(1,y);
			if(mode == 1) de[i] = -de[ip];
			else de[i] = de[ip];
			int j = cell(rowSize-1,y);
			int jp = cell(rowSize-2,y);
			if(mode == 1) de[j] = -de[jp];
			else de[j] = de[jp];
		}
		
		de[cell(0,0)] = 0.5f*(de[cell(1,0)]+de[cell(0,1)]);
		de[cell(0,colSize-1)] = 0.5f*(de[cell(1,colSize-1)]+de[cell(0,colSize-2)]);
		de[cell(rowSize-1,0)] = 0.5f*(de[cell(rowSize-2,0)]+de[cell(rowSize-1,1)]);
		de[cell(rowSize-1,colSize-1)] = 0.5f*(de[cell(rowSize-2,colSize-1)]+de[cell(rowSize-1,colSize-2)]);
		
		/***********************
		 * velocities for walls
		 ***********************/
		
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) {
					int pxu = ncell(x, y, x+1,y);
					int pxd = ncell(x, y, x-1,y);
					int pyu = ncell(x, y, x,y+1);
					int pyd = ncell(x, y, x,y-1);
					if(mode == 0) {
						float acc = 0;
						float div = 0;
						if(pxu != p) {
							acc += de[pxu];
							div += 1.0;
						}
						if(pxd != p) {
							acc += de[pxd];
							div += 1.0;
						}
						if(pyu != p) {
							acc += de[pyu];
							div += 1.0;
						}
						if(pyd != p) {
							acc += de[pyd];
							div += 1.0;
						}
						if(div > 0) de[p] = acc / div;
						else de[p] = 0;
					}
					if(mode == 1) {
						float acc = 0;
						float div = 0;
						if(pxu != p) {
							acc -= de[pxu];
							div += 1.0;
						}
						if(pxd != p) {
							acc -= de[pxd];
							div += 1.0;
						}
						if(pyu != p) {
							acc += de[pyu];
							div += 1.0;
						}
						if(pyd != p) {
							acc += de[pyd];
							div += 1.0;
						}
						if(div > 0) de[p] = acc / div;
						else de[p] = 0;
					}
					if(mode == 2) {
						float acc = 0;
						float div = 0;
						if(pxu != p) {
							acc += de[pxu];
							div += 1.0;
						}
						if(pxd != p) {
							acc += de[pxd];
							div += 1.0;
						}
						if(pyu != p) {
							acc -= de[pyu];
							div += 1.0;
						}
						if(pyd != p) {
							acc -= de[pyd];
							div += 1.0;
						}
						if(div > 0) de[p] = acc / div;
						else de[p] = 0;
					}
				}
			}
		}
	}
	
	public int cell(int x, int y) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		return x + y * rowSize;
	}
	
	public int ncell(int x, int y, int xn, int yn) {
		if(xn < 0) xn = 0;
		if(xn >= rowSize) xn = rowSize - 1;
		if(yn < 0) yn = 0;
		if(yn >= colSize) yn = colSize - 1;
		int pos = xn + yn * rowSize;
		if(wall[pos]) return cell(x, y);
		else return pos;
	}
}
