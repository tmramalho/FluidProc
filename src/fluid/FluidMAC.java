package fluid;

import processing.core.*;

public class FluidMAC {
	PApplet parent;
	private int cellw;
	private int cellh;
	private int rowSize;
	private int colSize;
	private int nw;
	private int nh;
	private int tsize;
	private int mxsize;
	private int mysize;
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
	private int[] reqCell;
	private boolean[] wall;
	private int[] numParticles;
	private float addVel;
	private boolean pv = true;
	private boolean dif = false;
	private float nt = 0; 
	private float baseDens;
	
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
	
	public FluidMAC(int cw, int ch, int w, int h, PApplet p, float bd) {
		rowSize = cw + 2;
		colSize = ch + 2;
		parent = p;
		cellw = Math.round(w / cw);
		cellh = Math.round(h / ch);
		nw = cw;
		nh = ch;
		
		tsize = rowSize*colSize;
		mxsize = (rowSize+1)*colSize;
		mysize = rowSize*(colSize+1);
		dens0 = new float[tsize];
		dens1 = new float[tsize];
		vx0 = new float[mxsize];
		vx1 = new float[mxsize];
		vx2 = new float[mxsize];
		vxs = new float[mxsize];
		vy0 = new float[mysize];
		vy1 = new float[mysize];
		vy2 = new float[mysize];
		vys = new float[mysize];
		req = new float[tsize];
		temp0 = new float[tsize];
		temp1 = new float[tsize];
		ink0 = new float[tsize];
		ink1 = new float[tsize];
		noise = new float[tsize];
		nx = new float[tsize];
		ny = new float[tsize];
		indReq = new float[4*tsize];
		reqCell = new int[4*tsize];
		wall = new boolean[tsize];
		numParticles = new int[tsize];
		addVel = 10;
		baseDens = bd;
		
		recalcNoiseFields(noise, nx, ny);
	}
	
	public void update(float a) {
		clearWall();
		
		wipe(vxs); wipe(vys);
		//pressure(dens1, vxs, vys, (float) (0.1*a));
		temp(temp1, vxs, vys, 0.001f*a);
		fall(ink1, vxs, vys, 0.0001f*a);
		addFields(vx0, vy0, vxs, vys);
		
		if(dif) diffuse(1, vx0, vx1, 0.0001f*a);
		else copy(vx0,vx1);
		if(dif)diffuse(2, vy0, vy1, 0.0001f*a);
		else copy(vy0,vy1);
		
		project(vx1, vy1, vx2, vy2);
		vconf(vx1, vy1, vx0, vy0, vx2, vy2, a);
		addTurbulence(vx1, vy1, nx, ny, a, 0.1f);
		setBnd(vx1, 1); setBnd(vy1, 2);
		
		copy(vx1,vx0);
		copy(vy1,vy0);
		copy(vx1,vx2);
		copy(vy1,vy2);
		
		sladvect(1, vx0, vy0, vx0, vx2, a);
		sladvect(2, vx0, vy0, vy0, vy2, a);
		
		//project(vx2, vy2, vx1, vy1);
		
		setBnd(vx2, 1); setBnd(vy2, 2);		
		
		copy(vx2, vx0);
		copy(vy2, vy0);
		
		if(dif) diffuse(0, dens0, dens1, 0.0001f*a);
		else copy(dens0, dens1);
		if(dif) diffuse(0, temp0, temp1, 0.001f*a);
		else copy(temp0, temp1);
		if(dif) diffuse(0, ink0, ink1, 0.0001f*a);
		else copy(ink0, ink1);
		
		sladvect(0, vx0, vy0, dens1, dens0, a); setBnd(dens0, 0); copy(dens0, dens1);
		sladvect(0, vx0, vy0, temp1, temp0, a); setBnd(temp0, 0); copy(temp0, temp1);
		sladvect(0, vx0, vy0, ink1, ink0, a); setBnd(ink0, 0); copy(ink0, ink1);
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
				parent.rect(x*cw, y*ch, cw, ch);
				if(pv) {
					float vxv = getMACCellX(x,y,vx0);
					float vyv = getMACCellY(x,y,vy0);
					if(Math.abs(vxv) > 0.01 || Math.abs(vyv) > 0.01) {
						parent.stroke(192, 0, 0);
						parent.line((float)(x*cw+0.5*cw), 
								(float)(y*ch+0.5*ch), 
								(float)(x*cw+0.5*cw+5*vxv), 
								(float)(y*ch+0.5*ch+5*vyv));
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
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		PApplet.println(getMACCellX(posx, posy, vx0) + " " 
				+ getMACCellY(posx, posy, vy0) + " " + getMACCellX(posx, posy, vx1) + " " 
				+ getMACCellY(posx, posy, vy1) + " " + getMACCellX(posx, posy, vx2) + " " 
				+ getMACCellY(posx, posy, vy2));
	}
	
	public void addDens(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float dens = dens0[p];
		dens = baseDens;
		dens0[p] = dens;
	}
	
	public void addTemp(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = temp0[p];
		temp = baseDens;
		temp0[p] = temp;
	}
	
	public void addInk(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float ink = ink0[p];
		ink = baseDens;
		ink0[p] = ink;
	}
	
	public void addDens(int mouseX, int mouseY, float bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float dens = dens0[p];
		dens += 10*bd;
		dens0[p] = dens;
	}
	
	public void addTemp(int mouseX, int mouseY, float bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = temp0[p];
		temp += bd;
		temp0[p] = temp;
	}
	
	public void addInk(int mouseX, int mouseY, float bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = ink0[p];
		temp += bd;
		ink0[p] = temp;
	}
	
	public void addVel(int mouseX, int mouseY, float angle) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		float nextX = (float) (addVel*Math.cos(angle));
		float nextY = (float) (addVel*Math.sin(angle));
		addtoMACCellX(posx, posy, vx0, nextX);
		addtoMACCellY(posx, posy, vy0, nextY);
	}
	
	public void addWall(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		wall[cell(posx, posy)] = true;
		if(posx == 1) wall[cell(0, posy)] = true;
		if(posx == rowSize-2) wall[cell(rowSize-1, posy)] = true;
		if(posy == 1) wall[cell(posx, 0)] = true;
		if(posy == colSize-2) wall[cell(posx, colSize-1)] = true;
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
	
	public void wipe(float[] d) {
		for(int i = 0; i < tsize; i++) d[i] = 0;
	}
	
	private void addFields(float[] dx, float[] dy, float[] sx, float[] sy) {
		for(int i = 0; i < mxsize; i++) {
			dx[i] += sx[i];
		}
		for(int i = 0; i < mysize; i++) {
			dy[i] += sy[i];
		}
	}
	
	public void reset() {
		for(int i = 0; i < tsize; i++) {
			wall[i] = false;
			dens0[i] = 0;
			dens1[i] = 0;
			temp0[i] = 0;
			temp1[i] = 0;
			ink0[i] = 0;
			ink1[i] = 0;
		}
		for(int i = 0; i < mxsize; i++) {
			vx0[i] = 0;
			vx1[i] = 0;
			vx2[i] = 0;
		}
		for(int i = 0; i < mysize; i++) {
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
				ax[pn] = (float) (0.5*nw*(an[pyu] - an[pyd]));
				ay[pn] = (float) (-0.5*nh*(an[pxu] - an[pxd]));
			}
		}
		
		nt += 10*nscale;
	}
	
	public void addTurbulence(float[] vx, float[] vy, float[] nx, float[] ny, float c, float d) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int i = cell(x,y);
				float vxi = getMACCellX(x, y, vx);
				float vyi = getMACCellY(x, y, vy);
				float len = (float) Math.sqrt(vxi*vxi + vyi*vyi);
				float scalex = (float) - x * (x - rowSize) / (float) (rowSize*rowSize);
				float scaley = (float) - y * (y - colSize) / (float) (colSize*colSize);
				if (len == 0) continue;
				float sc = len*len*len / (len*len*len + d);
				addtoMACCellX(x, y, vx, c*sc*scalex*scalex*nx[i]/nw);
				addtoMACCellY(x, y, vy, c*sc*scaley*scaley*ny[i]/nh);
			}
		}
		
		//recalcNoiseFields(noise, nx, ny);
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
	
	public void temp(float[] dens, float[] vx, float[] vy, float a) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				float ay = (float) Math.abs(dens[p]);
				addtoMACCellY(x,y,vy,-1.0f*a*ay);
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
	
	public void vconf(float[] vxa, float[] vya, float[] wx, float[] wy, float[] wa, float[] wu, float c) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				wx[p] = (float) (0.5*nw*(vya[cellY(x+1,y)] - vya[cellY(x,y)]));
				wy[p] = (float) (-0.5*nh*(vxa[cellX(x,y+1)] - vxa[cellX(x,y)]));
				wa[p] = (float) Math.sqrt(wx[p]*wx[p] + wy[p]*wy[p]);
			}
		}
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				float dwx = (float) (0.5*nw*(wa[pxu] - wa[pxd]));
				float dwy = (float) (0.5*nh*(wa[pyu] - wa[pyd]));
				float dn = 0;
				float n = (float) Math.sqrt(dwx*dwx+dwy*dwy);
				if(n != 0) {
					dn = (float) (1.0/n);
				}
				addtoMACCellX(x, y, vxa, c*(dwy * dn * wx[p])/nw);
				addtoMACCellY(x, y, vya, -1.0f*c*(dwx * dn * wy[p])/nh);
			}
		}
		
	}
	
	public void diffuse(int mode, float[] ain, float[] aout, float d) {
		float k = d*nh*nw;
		if(mode == 1) {
			for(int n=0; n< 5; n++) {
				wipe(aout);
				for(int x = 1; x < rowSize-1; x++) {
					for(int y = 1; y < colSize-1; y++) {
						int i = cell(x,y);
						if(wall[i]) continue;
						
						float res = (getMACCellX(x,y,ain) + k*(getMACCellX(x+1,y,aout) 
								+ getMACCellX(x-1,y,aout) 
								+ getMACCellX(x,y+1,aout) 
								+ getMACCellX(x,y-1,aout)))/(1+4*k);
						addtoMACCellX(x,y,aout, res);
					}
				}
				setBnd(aout, mode);
			}
		} else if(mode == 2) {
			for(int n=0; n< 5; n++) {
				wipe(aout);
				for(int x = 1; x < rowSize-1; x++) {
					for(int y = 1; y < colSize-1; y++) {
						int i = cell(x,y);
						if(wall[i]) continue;
						
						float res = (getMACCellY(x,y,ain) + k*(getMACCellY(x+1,y,aout) 
								+ getMACCellY(x-1,y,aout) 
								+ getMACCellY(x,y+1,aout) 
								+ getMACCellY(x,y-1,aout)))/(1+4*k);
						addtoMACCellY(x,y,aout, res);
					}
				}
				setBnd(aout, mode);
			}
		} else {
			for(int n=0; n< 5; n++) {
				for(int x = 1; x < rowSize-1; x++) {
					for(int y = 1; y < colSize-1; y++) {
						int i = cell(x,y);
						if(wall[i]) continue;
						int pxu = ncell(x, y, x+1,y);
						int pxd = ncell(x, y, x-1,y);
						int pyu = ncell(x, y, x,y+1);
						int pyd = ncell(x, y, x,y-1);
						aout[i] = (ain[i] + k*(aout[pxu] + aout[pxd] + aout[pyu] + aout[pyd]))/(1+4*k);
					}
				}
				setBnd(aout, mode);
			}
		}
	}
	
	public void project(float[] vx, float[] vy, float[] div, float[] p) {
		for(int x = 1; x < rowSize-1; x++) {
			for(int y = 1; y < colSize-1; y++) {
				int i = cell(x,y);
				div[i] = -0.5f*((vx[cellX(x+1,y)] - vx[cellX(x,y)])/nw + (vy[cellY(x, y+1)] - vy[cellY(x, y)])/nh);
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
				int pxu = cell(x+1, y);
				int pxd = cell(x-1, y);
				int pyu = cell(x, y+1);
				int pyd = cell(x, y-1);
				addtoMACCellX(x, y, vx, -0.5f*nw*(p[pxu]-p[pxd]));
				addtoMACCellY(x, y, vy, -0.5f*nh*(p[pyu]-p[pyd]));
			}
		}	
	}
	
	public void sladvect(int mode, float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		wipe(aout);
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				
				float vx = getMACCellX(x, y, vxa);
				float vy = getMACCellY(x, y, vya);
					
				float xn = x - vx * c * nw;
				float yn = y - vy * c * nh;
				
				float xb = (float) (rowSize - 1 + 0.5);
				float yb = (float) (colSize - 1 + 0.5);
				
				if(xn < 0.5) xn = (float) 0.5;
				else if (xn > xb) xn = xb;
				if(yn < 0.5) yn = (float) 0.5;
				else if (yn > yb) yn = yb;
				
				int x0 = (int) xn;
				int y0 = (int) yn;
				int x1 = x0 + 1;
				int y1 = y0 + 1;
				
				float s1 = xn - x0;
				float s0 = 1 - s1;
				float t1 = yn - y0;
				float t0 = 1 - t1;
				
				if(mode == 0) {
					aout[p] = s0*(t0*ain[cell(x0, y0)] + t1*ain[cell(x0, y1)]) + 
							  s1*(t0*ain[cell(x1, y0)] + t1*ain[cell(x1, y1)]);
				} else if (mode == 1) {
					float res = s0*(t0*getMACCellX(x0, y0, ain) + t1*getMACCellX(x0, y1, ain)) + 
					  			s1*(t0*getMACCellX(x1, y0, ain) + t1*getMACCellX(x1, y1, ain));
					addtoMACCellX(x, y, aout, res);
				} else if (mode == 2) {
					float res = s0*(t0*getMACCellY(x0, y0, ain) + t1*getMACCellY(x0, y1, ain)) + 
		  						s1*(t0*getMACCellY(x1, y0, ain) + t1*getMACCellY(x1, y1, ain));
					addtoMACCellY(x, y, aout, res);
				}
			}
		}
	}
	
	private void setBnd(float[] de, int mode) {
		if(mode == 1) {
			for(int x = 1; x < rowSize; x++) {
				de[cellX(x,0)] = de[cellX(x,1)];
				de[cellX(x,colSize-1)] = de[cellX(x,colSize-2)];
			}
			for(int y = 1; y < colSize - 1; y++) {
				de[cellX(0,y)] = -de[cellX(1,y)];
				de[cellX(rowSize,y)] = -de[cellX(rowSize-1,y)];
			}
			de[cellX(0,0)] = 0.5f*(de[cellX(1,0)]+de[cellX(0,1)]);
			de[cellX(0,colSize-1)] = 0.5f*(de[cellX(1,colSize-1)]+de[cellX(0,colSize-2)]);
			de[cellX(rowSize,0)] = 0.5f*(de[cellX(rowSize-1,0)]+de[cellX(rowSize,1)]);
			de[cellX(rowSize,colSize-1)] = 0.5f*(de[cellX(rowSize-1,colSize-1)]+de[cellX(rowSize,colSize-2)]);
		} else if(mode == 2) {
			for(int x = 1; x < rowSize - 1; x++) {
				de[cellY(x,0)] = -de[cellY(x,1)];
				de[cellY(x,colSize)] = -de[cellY(x,colSize-1)];
			}
			for(int y = 1; y < colSize; y++) {
				de[cellY(0,y)] = de[cellY(1,y)];
				de[cellY(rowSize-1,y)] = de[cellY(rowSize-2,y)];
			}
			de[cellY(0,0)] = 0.5f*(de[cellY(1,0)]+de[cellY(0,1)]);
			de[cellY(0,colSize)] = 0.5f*(de[cellY(1,colSize)]+de[cellY(0,colSize-1)]);
			de[cellY(rowSize-1,0)] = 0.5f*(de[cellY(rowSize-2,0)]+de[cellY(rowSize-1,1)]);
			de[cellY(rowSize-1,colSize)] = 0.5f*(de[cellY(rowSize-2,colSize)]+de[cellY(rowSize-1,colSize-1)]);
		} else {
			for(int x = 1; x < rowSize - 1; x++) {
				de[cell(x,0)] = de[cell(x,1)];
				de[cell(x,colSize-1)] = de[cell(x,colSize-2)];
			}
			for(int y = 1; y < colSize - 1; y++) {
				de[cell(0,y)] = de[cell(1,y)];
				de[cell(rowSize-1,y)] = de[cell(rowSize-2,y)];
			}
			de[cell(0,0)] = 0.5f*(de[cell(1,0)]+de[cell(0,1)]);
			de[cell(0,colSize-1)] = 0.5f*(de[cell(1,colSize-1)]+de[cell(0,colSize-2)]);
			de[cell(rowSize-1,0)] = 0.5f*(de[cell(rowSize-2,0)]+de[cell(rowSize-1,1)]);
			de[cell(rowSize-1,colSize-1)] = 0.5f*(de[cell(rowSize-2,colSize-1)]+de[cell(rowSize-1,colSize-2)]);
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
						if(pxu != p) {
							de[cellX(x+1,y)] = -de[cellX(x+2,y)];
						}
						if(pxd != p) {
							de[cellX(x,y)] = -de[cellX(x-1,y)];
						}
					}
					if(mode == 2) {
						if(pyu != p) {
							de[cellY(x,y+1)] = -de[cellY(x,y+2)];
						}
						if(pyd != p) {
							de[cellY(x,y)] = -de[cellY(x,y-1)];
						}
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
	
	public int cellX(int x, int y) {
		return x + y * (rowSize + 1);
	}
	
	public int cellY(int x, int y) {
		return x + y * rowSize;
	}
	
	public void addtoMACCellX(int x, int y, float[] arr, float val) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellX(x , y);
		int pn = cellX(x+1 , y);
		arr[p] += val * 0.5f;
		arr[pn] += val * 0.5f;
	}
	
	public void addtoMACCellY(int x, int y, float[] arr, float val) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellY(x , y);
		int pn = cellY(x , y+1);
		arr[p] += val * 0.5f;
		arr[pn] += val * 0.5f;
	}
	
	public float getMACCellX(int x, int y, float[] arr) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellX(x , y);
		int pn = cellX(x+1 , y);
		return (arr[p]+arr[pn]) * 0.5f;
	}
	
	public float getMACCellY(int x, int y, float[] arr) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellY(x , y);
		int pn = cellY(x , y+1);
		return (arr[p]+arr[pn]) * 0.5f;
	}
	
	/*****************************
	 * Code for compressible flow
	 *****************************/
	
	public void cproject(float[] vx, float[] vy, float[] div, float[] p, float[] sx, float[] sy, float c) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int i = cell(x,y);
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				div[i] = (float) (-0.5*((vx[pxu] - vx[pxd])/nw + (vy[pyu] - vy[pyd])/nh));
				div[i] += (float) c*(-0.5*((sx[pxu] - sx[pxd])/nw + (sy[pyu] - sy[pyd])/nh));
				p[i] = 0;
			}
		}
		setBnd(p, 0);
		setBnd(div, 0);
		
		for(int n=0; n< 20; n++) {
			for(int x = 1; x < rowSize-1; x++) {
				for(int y = 1; y < colSize-1; y++) {
					int i = cell(x,y);
					int pxu = ncell(x, y, x+1,y);
					int pxd = ncell(x, y, x-1,y);
					int pyu = ncell(x, y, x,y+1);
					int pyd = ncell(x, y, x,y-1);
					p[i] = (div[i] + p[pxu] + p[pxd] + p[pyu] + p[pyd])/4;
				}
			}
			setBnd(p, 0);
		}
		
		for(int x = 1; x < rowSize-1; x++) {
			for(int y = 1; y < colSize-1; y++) {
				int i = cell(x,y);
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				vx[i] -= 0.5*nw*(p[pxu]-p[pxd]);
				vy[i] -= 0.5*nh*(p[pyu]-p[pyd]);
			}
		}	
	}
	
	public void sluadvect(float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		wipe(req);
		
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				
				float vx = vxa[p];
				float vy = vya[p];
					
				float xn = x - vx * c * nw;
				float yn = y - vy * c * nh;
				
				float xb = (float) (rowSize - 1 + 0.5);
				float yb = (float) (colSize - 1 + 0.5);
				
				if(xn < 0.5) xn = (float) 0.5;
				else if (xn > xb) xn = xb;
				if(yn < 0.5) yn = (float) 0.5;
				else if (yn > yb) yn = yb;
				
				int x0 = (int) xn;
				int y0 = (int) yn;
				int x1 = x0 + 1;
				int y1 = y0 + 1;
				
				float s1 = xn - x0;
				float s0 = 1 - s1;
				float t1 = yn - y0;
				float t0 = 1 - t1;
				
				float ia = s0*t0;
				float ib = s0*t1;
				float ic = s1*t0;
				float id = s1*t1;
				
				//save sources and what each loses
				reqCell[4*p]   = cell(x0, y0);
				reqCell[4*p+1] = cell(x0, y1);
				reqCell[4*p+2] = cell(x1, y0);
				reqCell[4*p+3] = cell(x1, y1);
				
				indReq[4*p] = ia;
				indReq[4*p+1] = ib;
				indReq[4*p+2] = ic;
				indReq[4*p+3] = id;
				
				//accumulate how much each cell loses in total
				req[cell(x0, y0)] += ia;
				req[cell(x0, y1)] += ib;
				req[cell(x1, y0)] += ic;
				req[cell(x1, y1)] += id;
			}
		}
		
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				//recover previous data
				float ia = indReq[4*p];
				float ib = indReq[4*p+1];
				float ic = indReq[4*p+2];
				float id = indReq[4*p+3];
				
				//recover source cells
				int sc1 = reqCell[4*p];
				int sc2 = reqCell[4*p+1];
				int sc3 = reqCell[4*p+2];
				int sc4 = reqCell[4*p+3];
				
				//get total fractions and rescale requests
				float fa = req[sc1];
				float fb = req[sc2];
				float fc = req[sc3];
				float fd = req[sc4];
				
				if (fa<1.0f) fa = 1.0f;
				if (fb<1.0f) fb = 1.0f;
				if (fc<1.0f) fc = 1.0f;
				if (fd<1.0f) fd = 1.0f;
				
				ia = ia * ain[sc1] / fa;
				ib = ib * ain[sc2] / fb;
				ic = ic * ain[sc3] / fc;
				id = id * ain[sc4] / fd;
				
				aout[p] += (ia+ib+ic+id);
				aout[sc1] -= ia;
				aout[sc2] -= ib;
				aout[sc3] -= ic;
				aout[sc4] -= id;		
			}
		}
	}
}
