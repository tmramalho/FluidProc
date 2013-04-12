package fluid;

import processing.core.*;

public class FluidCom {
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
	private float[] fuel0;
	private float[] fuel1;
	private float[] ink0;
	private float[] ink1;
	private float[] indReq;
	private float[] noise;
	private float[] nx;
	private float[] ny;
	private int[] reqCell;
	private int[] lastCoord;
	private boolean[] wall;
	private int[] numParticles;
	private boolean pv = false;
	private boolean dif = false;
	private float nt = 0; 
	private float vol = 0.1f; //volatility
	private float addVel = 100;
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
	
	public FluidCom(int cw, int ch, int w, int h, PApplet p, float bd) {
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
		fuel0 = new float[tsize];
		fuel1 = new float[tsize];
		ink0 = new float[tsize];
		ink1 = new float[tsize];
		noise = new float[tsize];
		nx = new float[tsize];
		ny = new float[tsize];
		indReq = new float[4*tsize];
		reqCell = new int[4*tsize];
		lastCoord = new int[2];
		wall = new boolean[tsize];
		numParticles = new int[tsize];
		baseDens = bd;
		
		recalcNoiseFields(noise, nx, ny);
	}
	
	public void update(float a) {
		clearWall();
		
		wipe(vxs); wipe(vys);
		//pressure(dens1, vxs, vys, 0.1f*a);
		temp(temp1, vxs, vys, 0.01f*a);
		fall(ink1, vxs, vys, 0.0001f*a);
		addFields(vx0, vy0, vxs, vys);
		
		if(dif) diffuse(1, vx0, vx1, 0.0001f*a);
		else copy(vx0,vx1);
		if(dif)diffuse(2, vy0, vy1, 0.0001f*a);
		else copy(vy0,vy1);
		
		//project(vx1, vy1, vx2, vy2);
		cproject(vx1, vy1, vx2, vy2, fuel0, vol);
		vconf(vx1, vy1, vx0, vy0, vx2, vy2, 100*a);
		//addTurbulence(vx1, vy1, nx, ny, a, 0.1f);
		setBnd(vx1, 1); setBnd(vy1, 2);
		
		copy(vx1,vx0);
		copy(vy1,vy0);
		copy(vx1,vx2);
		copy(vy1,vy2);
		//copy(vx0, vxs); //debug only
		//copy(vy0, vys); //debug only
		
		sladvect(1, vx0, vy0, vx0, vx2, a);
		sladvect(2, vx0, vy0, vy0, vy2, a);
		revsladvect(1, vx2, vy2, vx2, vx1, a);
		revsladvect(2, vx2, vy2, vy2, vy1, a);
		mcCorrect(1, vx0, vy0, vx1, vx2, a);
		mcCorrect(2, vx0, vy0, vy1, vy2, a);
		
		//project(vx2, vy2, vx1, vy1);
		cproject(vx2, vy2, vx0, vy0, fuel0, vol);
		
		setBnd(vx2, 1); setBnd(vy2, 2);		
		
		copy(vx2, vx0);
		copy(vy2, vy0);
		
		if(dif) { diffuse(0, dens0, dens1, 0.0001f*a); copy(dens1, dens0); }
		else copy(dens0, dens1);
		sluadvect(vx0, vy0, dens1, dens0, a); /*setBnd(dens0, 0);*/ copy(dens0, dens1);
		
		if(dif) { diffuse(0, temp0, temp1, 0.001f*a); copy(temp1, temp0); }
		else copy(temp0, temp1);
		sluadvect(vx0, vy0, temp1, temp0, a); /*setBnd(temp0, 0);*/ copy(temp0, temp1);
		
		if(dif) { diffuse(0, fuel0, fuel1, 0.001f*a); copy(fuel1, fuel0); }
		else copy(fuel0, fuel1);
		sluadvect(vx0, vy0, fuel1, fuel0, a); /*setBnd(ink0, 0);*/ copy(fuel0, fuel1);
		
		if(dif) { diffuse(0, ink0, ink1, 0.0001f*a); copy(ink1, ink0); }
		else copy(ink0, ink1);
		sluadvect(vx0, vy0, ink1, ink0, a); /*setBnd(ink0, 0);*/ copy(ink0, ink1);
		
		coolDown(fuel0, a);
		coolDown(temp0, a);
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
			}
		}
		if(pv) {
			for(int x = 1; x < rowSize - 1; x++) {
				for(int y = 1; y < colSize - 1; y++) {
					float cx = (x-1)*cellw+0.5f*cellw;
					float cy = (y-1)*cellh+0.5f*cellh;
					float vx = getMACCellX(x,y,vx0);
					float vy = getMACCellY(x,y,vy0);
					if(Math.abs(vx) > 0.01 || Math.abs(vy) > 0.01) {
						parent.stroke(0,192,0);
						parent.line(cx,cy,cx+3*vx,cy+3*vy);
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
				int fuel = (int) (0.5*fuel0[p]);
				int wallc = 0;
				if(wall[p]) wallc = 255;
				parent.fill(dens + temp + fuel, dens + wallc + fuel, dens + ink);
				parent.rect(x*cw, y*ch, cw, ch);
			}
		}
		if(pv) {
			for(int x = 0; x < rowSize; x++) {
				for(int y = 0; y < colSize; y++) {
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
					float vxvd = getMACCellX(x,y,vx1);
					float vyvd = getMACCellY(x,y,vy1);
					if(Math.abs(vxvd) > 0.01 || Math.abs(vyvd) > 0.01) {
						parent.stroke(0, 128, 0);
						parent.line((float)(x*cw+0.5*cw), 
								(float)(y*ch+0.5*ch), 
								(float)(x*cw+0.5*cw+5*vxvd), 
								(float)(y*ch+0.5*ch+5*vyvd));
						parent.stroke(0);
					}
					float vxvs = getMACCellX(x,y,vxs);
					float vyvs = getMACCellY(x,y,vys);
					if(Math.abs(vxvs) > 0.01 || Math.abs(vyvs) > 0.01) {
						parent.stroke(128, 128, 0);
						parent.line((float)(x*cw+0.5*cw), 
								(float)(y*ch+0.5*ch), 
								(float)(x*cw+0.5*cw+5*vxvs), 
								(float)(y*ch+0.5*ch+5*vyvs));
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
		PApplet.println(getMACCellX(posx, posy, vx0) + ": " 
				+ vx0[cellX(posx, posy)] + " " + vx0[cellX(posx + 1, posy)] + " ; " 
				+ getMACCellY(posx, posy, vy0) + ": " 
				+ vy0[cellY(posx, posy)] + " " + vy0[cellY(posx, posy + 1)]);
	}
	
	public void printAccs(String s) {
		float tacc = 0;
		float iacc = 0;
		float dacc = 0;
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				tacc += temp0[cell(x,y)];
				iacc += ink0[cell(x,y)];
				dacc += dens0[cell(x,y)];
			}
		}
		PApplet.println(s + " Temp " + tacc + " Dens " + dacc + " Ink " + iacc);
	}
	
	public void addDens(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float dens = dens0[p];
		dens = baseDens;
		dens0[p] = dens;
	}
	
	public void addLinDens(int mouseXf, int mouseYf, int mouseX, int mouseY) {
		int x1 = (int)(mouseX / cellw) + 1;
		int y1 = (int)(mouseY / cellh) + 1;
		int x0 = (int)(mouseXf / cellw) + 1;
		int y0 = (int)(mouseYf / cellh) + 1;
				
		int dx = Math.abs(x1-x0);
		int dy = Math.abs(y1-y0);
		int sx = -1;
		int sy = -1;
		if (x0 < x1) sx = 1;
		if (y0 < y1) sy = 1;
		int err = dx-dy;
		 
		while(true) {
			dens0[cell(x0,y0)] = baseDens;
			if (x0 == x1 && y0 == y1) break;
			int e2 = 2*err;
			if (e2 > -dy) { 
				err = err - dy;
				x0 = x0 + sx;
			}
			if (e2 < dx) { 
				err = err + dx;
				y0 = y0 + sy;
			}
		}
	}
	
	public void addTemp(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float temp = temp0[p];
		temp = baseDens;
		temp0[p] = temp;
	}
	
	public void addFuel(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		float fuel = fuel0[p];
		fuel = baseDens;
		fuel0[p] = fuel;
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
			if(posx == 1) wall[cell(0, posy)] = false;
			if(posx == rowSize-2) wall[cell(rowSize-1, posy)] = false;
			if(posy == 1) wall[cell(posx, 0)] = false;
			if(posy == colSize-2) wall[cell(posx, colSize-1)] = false;
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
			fuel0[i] = 0;
			fuel1[i] = 0;
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
				
				vy[cellY(x,y)] += a*(dens[p] - dens[pyd]);
				vy[cellY(x,y+1)] += a*(dens[pyu] - dens[p]);
				vx[cellX(x,y)] += a*(dens[p] - dens[pxd]);
				vx[cellX(x+1,y)] += a*(dens[pxu] - dens[p]);
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
				float dx = vx * c * nw;
				float dy = vy * c * nh;
				float xn = 0;
				float yn = 0;
					
				if(Math.abs(dx) > 1 || Math.abs(dy) > 1) {
					pathCollisions(x, y, (int) (x - dx), (int) (y - dy));
					xn = lastCoord[0];
					yn = lastCoord[1];
				} else {
					xn = x - dx;
					yn = y - dy;
				}
				
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
	
	public void revsladvect(int mode, float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		wipe(aout);
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				
				float vx = getMACCellX(x, y, vxa);
				float vy = getMACCellY(x, y, vya);
				float dx = vx * c * nw;
				float dy = vy * c * nh;
				float xn = 0;
				float yn = 0;
					
				if(Math.abs(dx) > 1 || Math.abs(dy) > 1) {
					pathCollisions(x, y, (int) (x + dx), (int) (y + dy));
					xn = lastCoord[0];
					yn = lastCoord[1];
				} else {
					xn = x + dx;
					yn = y + dy;
				}
				
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
	
	public void mcCorrect(int mode, float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				
				float vx = getMACCellX(x, y, vxa);
				float vy = getMACCellY(x, y, vya);
				float dx = vx * c * nw;
				float dy = vy * c * nh;
				float xn = 0;
				float yn = 0;
					
				if(Math.abs(dx) > 1 || Math.abs(dy) > 1) {
					pathCollisions(x, y, (int) (x - dx), (int) (y - dy));
					xn = lastCoord[0];
					yn = lastCoord[1];
				} else {
					xn = x - dx;
					yn = y - dy;
				}
				
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
					float dv = (getMACCellX(x, y, vxa) - getMACCellX(x, y, ain))/2;
					float res = Math.abs(getMACCellX(x, y, aout) + dv);
					
					float min = getMACCellX(x0, y0, vxa);
					if (getMACCellX(x0, y1, vxa) < min ) min = getMACCellX(x0, y1, vxa);
					if (getMACCellX(x1, y0, vxa) < min ) min = getMACCellX(x1, y0, vxa);
					if (getMACCellX(x1, y1, vxa) < min ) min = getMACCellX(x1, y1, vxa);
					
					float max = getMACCellX(x0, y0, vxa);
					if (getMACCellX(x0, y1, vxa) > max ) max = getMACCellX(x0, y1, vxa);
					if (getMACCellX(x1, y0, vxa) > max ) max = getMACCellX(x1, y0, vxa);
					if (getMACCellX(x1, y1, vxa) > max ) max = getMACCellX(x1, y1, vxa);
					
					if(res > Math.abs(min) && res < Math.abs(max)){
						addtoMACCellX(x, y, aout, dv);
					}
				} else if (mode == 2) {
					float dv = (getMACCellY(x, y, vya) - getMACCellY(x, y, ain))/2;
					float res = Math.abs(getMACCellY(x, y, aout) + dv);
					
					float min = getMACCellY(x0, y0, vya);
					if (getMACCellY(x0, y1, vya) < min ) min = getMACCellY(x0, y1, vya);
					if (getMACCellY(x1, y0, vya) < min ) min = getMACCellY(x1, y0, vya);
					if (getMACCellY(x1, y1, vya) < min ) min = getMACCellY(x1, y1, vya);
					
					float max = getMACCellY(x0, y0, vya);
					if (getMACCellY(x0, y1, vya) > max ) max = getMACCellY(x0, y1, vya);
					if (getMACCellY(x1, y0, vya) > max ) max = getMACCellY(x1, y0, vya);
					if (getMACCellY(x1, y1, vya) > max ) max = getMACCellY(x1, y1, vya);
					
					if(res > Math.abs(min) && res < Math.abs(max)){
						addtoMACCellY(x, y, aout, dv);
					}
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
						} else de[cellX(x+1,y)] = 0;
						if(pxd != p) {
							de[cellX(x,y)] = -de[cellX(x-1,y)];
						} else de[cellX(x,y)] = 0;
					}
					if(mode == 2) {
						if(pyu != p) {
							de[cellY(x,y+1)] = -de[cellY(x,y+2)];
						} else de[cellY(x,y+1)] = 0;
						if(pyd != p) {
							de[cellY(x,y)] = -de[cellY(x,y-1)];
						} else de[cellY(x,y)] = 0;
					}
				}
			}
		}
	}
	
	final public void imposeWall(float[] dp) {
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				if(wall[cell(x,y)])  dp[cell(x,y)] = 0;
			}
		}
	}
	
	final public int cell(int x, int y) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		return x + y * rowSize;
	}
	
	final public int ncell(int x, int y, int xn, int yn) {
		if(xn < 0) xn = 0;
		if(xn >= rowSize) xn = rowSize - 1;
		if(yn < 0) yn = 0;
		if(yn >= colSize) yn = colSize - 1;
		int pos = xn + yn * rowSize;
		if(wall[pos]) return cell(x, y);
		else return pos;
	}
	
	final public int cellX(int x, int y) {
		return x + y * (rowSize + 1);
	}
	
	final public int cellY(int x, int y) {
		return x + y * rowSize;
	}
	
	final public void addtoMACCellX(int x, int y, float[] arr, float val) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellX(x , y);
		int pn = cellX(x+1 , y);
		arr[p] += val * 0.5f;
		arr[pn] += val * 0.5f;
	}
	
	final public void addtoMACCellY(int x, int y, float[] arr, float val) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellY(x , y);
		int pn = cellY(x , y+1);
		arr[p] += val * 0.5f;
		arr[pn] += val * 0.5f;
	}
	
	final public float getMACCellX(int x, int y, float[] arr) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		int p = cellX(x , y);
		int pn = cellX(x+1 , y);
		return (arr[p]+arr[pn]) * 0.5f;
	}
	
	final public float getMACCellY(int x, int y, float[] arr) {
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
	
	public void cproject(float[] vx, float[] vy, float[] div, float[] p, float[] tmp, float v) {
		for(int x = 1; x < rowSize-1; x++) {
			for(int y = 1; y < colSize-1; y++) {
				int i = cell(x,y);
				div[i] = -0.5f*((vx[cellX(x+1,y)] - vx[cellX(x,y)])/nw + (vy[cellY(x, y+1)] - vy[cellY(x, y)])/nh);
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				float dt = 0.5f*((tmp[pxu] - tmp[pxd])/nw + (tmp[pyu] - tmp[pyd])/nh);
				dt = dt < 1 ? 0 : dt;
				div[i] += v*dt;
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
			imposeWall(p);
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
	
	public void coolDown(float[] res, float dt) {
		float tacc = 0;
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				tacc += res[cell(x,y)];
			}
		}
		tacc *= 1.0f/(nw*nh);

		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				if(res[cell(x,y)] - tacc*dt < 0) {
					res[cell(x,y)] = 0;
				} else {
					res[cell(x,y)] -= tacc*dt;
				}
			}
		}
	}
	
	public void sluadvect(float[] vxa, float[] vya, float[] ain, float[] aout, float c) {
		wipe(req);
		
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				
				float vx = getMACCellX(x, y, vxa);
				float vy = getMACCellY(x, y, vya);
				float dx = vx * c * nw;
				float dy = vy * c * nh;
				float xn = 0;
				float yn = 0;
				
				if(Math.abs(dx) > 1 || Math.abs(dy) > 1) {
					pathCollisions(x, y, (int) (x - dx), (int) (y - dy));
					xn = lastCoord[0];
					yn = lastCoord[1];
				} else {
					xn = x - dx;
					yn = y - dy;
				}
				
				float xb = (float) (rowSize - 1 + 0.5);
				float yb = (float) (colSize - 1 + 0.5);
				
				if(xn < 0.5f) xn = 0.5f;
				else if (xn > xb) xn = xb;
				if(yn < 0.5f) yn = 0.5f;
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
	
	private void pathCollisions(int x0, int y0, int x1, int y1) {
		int dx = Math.abs(x1-x0);
		int dy = Math.abs(y1-y0);
		int sx = -1;
		int sy = -1;
		if (x0 < x1) sx = 1;
		if (y0 < y1) sy = 1;
		int err = dx-dy;
		int x0p = x0;
		int y0p = y0;
		 
		while(true) {
			if (x0 == x1 && y0 == y1) {
				//finished the line without hitting something
				lastCoord[0] = x1;
				lastCoord[1] = y1;
				break;
			}
			int e2 = 2*err;
			if (e2 > -dy) { 
				err = err - dy;
				x0p = x0;
				x0 = x0 + sx;
			}
			if (e2 < dx) { 
				err = err + dx;
				y0p = y0;
				y0 = y0 + sy;
			}
			if(wall[cell(x0,y0)] || wall[cell(x0,y0)] ) {
				//something is in the way
				//so return previous open point
				lastCoord[0] = x0p;
				lastCoord[1] = y0p;
				break;
			}
		}
	}
}

/*fixing walls?
if(mode == 0) {
	int p1 = cell(x0, y0); if(wall[p1]) p1 = p;
	int p2 = cell(x0, y1); if(wall[p2]) p2 = p;
	int p3 = cell(x1, y0); if(wall[p3]) p3 = p;
	int p4 = cell(x1, y1); if(wall[p4]) p4 = p;
	aout[p] = s0*(t0*ain[p1] + t1*ain[p2]) + 
			  s1*(t0*ain[p3] + t1*ain[p4]);
} else if (mode == 1) {
	float f1 = getMACCellX(x0, y0, ain); int p1 = cell(x0, y0); if(wall[p1]) f1 = getMACCellX(x, y, ain);
	float f2 = getMACCellX(x0, y1, ain); int p2 = cell(x0, y1); if(wall[p2]) f2 = getMACCellX(x, y, ain);
	float f3 = getMACCellX(x1, y0, ain); int p3 = cell(x1, y0); if(wall[p3]) f3 = getMACCellX(x, y, ain);
	float f4 = getMACCellX(x1, y1, ain); int p4 = cell(x1, y1); if(wall[p4]) f4 = getMACCellX(x, y, ain);
	float res = s0*(t0*f1 + t1*f2) + 
	  			s1*(t0*f3 + t1*f4);
	addtoMACCellX(x, y, aout, res);
} else if (mode == 2) {
	float f1 = getMACCellY(x0, y0, ain); int p1 = cell(x0, y0); if(wall[p1]) f1 = getMACCellY(x, y, ain);
	float f2 = getMACCellY(x0, y1, ain); int p2 = cell(x0, y1); if(wall[p2]) f2 = getMACCellY(x, y, ain);
	float f3 = getMACCellY(x1, y0, ain); int p3 = cell(x1, y0); if(wall[p3]) f3 = getMACCellY(x, y, ain);
	float f4 = getMACCellY(x1, y1, ain); int p4 = cell(x1, y1); if(wall[p4]) f4 = getMACCellY(x, y, ain);
	float res = s0*(t0*f1 + t1*f2) + 
	  			s1*(t0*f3 + t1*f4);
	addtoMACCellY(x, y, aout, res);
}*/
