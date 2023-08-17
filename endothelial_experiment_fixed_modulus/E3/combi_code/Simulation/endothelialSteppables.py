from cc3d.core.PySteppables import *
from fipy import CellVariable, Grid2D, dump
from test import tester
import variablevals as vv
import gc
import os
import sys
import random
import time
global mesh
from numpy import *
import numpy as np
from builtins import range
global cytokines
global cellpresent
global fullFileName
import csv

# globals for search
global setLambda
global setSaturationCoef
# global sigmoida
# global sigmoidb


# sigmoida = 
# sigmoidb = 


setlambda = 2000 #2000
setSaturationCoef = 10**-11
endocount = 500

def randopos(xboundary, yboundary, typein):
    if typein!=1:
        randxhi = random.randint(vv.nx-xboundary+2, vv.nx-2)
        randxlo = random.randint(2, xboundary-2)
        randyall = random.randint(2, vv.ny-2)
        randyhi = random.randint(vv.ny-yboundary+2, vv.ny-2)
        randylo = random.randint(2, yboundary-2)
        randxall = random.randint(2, vv.nx-2)
        hixvals = [randxhi, randyall]
        loxvals = [randxlo, randyall]
        hiyvals = [randxall, randyhi]
        loyvals = [randxall, randylo]
    if typein == 1:
        randxhi = random.randint(xboundary+2, vv.nx-2-xboundary)
        randxlo = random.randint(xboundary+2, vv.nx-2-xboundary)
        randyall = random.randint(yboundary+2, vv.ny-2-yboundary)
        randyhi = random.randint(yboundary+2, vv.ny-2-yboundary)
        randylo = random.randint(yboundary+2, vv.ny-2-yboundary)
        randxall = random.randint(xboundary+2, vv.nx-2-xboundary)
        hixvals = [randxhi, randyall]
        loxvals = [randxlo, randyall]
        hiyvals = [randxall, randyhi]
        loyvals = [randxall, randylo]
    rands = [hixvals, loxvals, hiyvals, loyvals]
    randselect = random.randint(0,4)
    randx = rands[randselect][0]
    randy = rands[randselect][1]
    # print(randx, randy, typein)
    return [randx, randy]
    
def sigmoid(x, a, b):
  
    z = np.exp(-a*(x - b))
    sig = 1 / (1 + z)

    return sig
    
def modmm(x, km, vmax):
    mm = abs((vmax * x)/(km + x))
    
    return mm

nx = vv.nx
ny = nx
dx = vv.dx
dy = dx
L = dx * nx
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
cellpresente = CellVariable(name = "solution variable",
           mesh = mesh,
           value = 0)
cellpresentndn = CellVariable(name = "solution variable",
           mesh = mesh,
           value = 0)
cellpresentna = CellVariable(name = "solution variable",
           mesh = mesh,
           value = 0)
cellpresentm1 = CellVariable(name = "solution variable",
           mesh = mesh,
           value = 0)
cellpresentm2 = CellVariable(name = "solution variable",
           mesh = mesh,
           value = 0)      
cytokines = [CellVariable(name = "solution variable", mesh = mesh, value = 0.0) for i in range(vv.total_cytokines)]

class endothelialSteppable(SteppableBasePy):

    def __init__(self,frequency=int(vv.relaxationmcs)):
        SteppableBasePy.__init__(self,frequency)
     

    def start(self):
        global fullFileName   
        
        global setSaturationCoef
        global endocount
        self.scalarFieldil8=self.create_scalar_field_py("il8")
        self.scalarFieldil1=self.create_scalar_field_py("il1")
        self.scalarFieldil6=self.create_scalar_field_py("il6")
        self.scalarFieldil10=self.create_scalar_field_py("il10")
        self.scalarFieldtnf=self.create_scalar_field_py("tnf")
        self.scalarFieldtgf=self.create_scalar_field_py("tgf")
        
        fileHandle,fullFileName=self.open_file_in_simulation_output_folder("datafiles/creatdoc.txt","w")
        if not os.path.exists(os.path.dirname(fullFileName)):
            os.makedirs(os.path.dirname(fullFileName))
        self.plot_win = self.add_new_plot_window(title='Cell counts',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell types', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        

        self.plot_win.add_plot('Endothelial', style='Lines',color='green', alpha=100, size=5)
        self.plot_win.add_plot('Neutrophils', style='Lines',color='red', alpha=100, size=5)
        self.plot_win.add_plot('Monocytes', style='Lines',color='blue', alpha=100, size=5)
        self.plot_win.add_plot('Fibroblast', style='Lines',color='orange', alpha=100, size=5)
        self.plot_win.add_plot('Neutrophil A', style='Lines',color='pink', alpha=100, size=5)
        self.plot_win.add_plot('Neutrophil NDN', style='Lines',color='white', alpha=100, size=5)
        self.plot_win.add_plot('Monocyte R', style='Lines',color='violet', alpha=100, size=5)
        self.plot_win.add_plot('Macrophage I', style='Lines',color='brown', alpha=100, size=5)
        self.plot_win.add_plot('Macrophage II', style='Lines',color='cyan', alpha=100, size=5)
        self.plot_win.add_plot('Myofibroblast', style='Lines',color='yellow', alpha=100, size=5)
        """
        any code in the start function runs before MCS=0
        """
        # size of cell will be 3x3x1
        
        for i in range(0,endocount):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.ENDOTHELIAL)
        for i in range(0,1000):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 2)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.NEUTROPHIL)
        for i in range(0,100):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.NEUTROPHIL)
        for i in range(0,100):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.FIBROBLAST)
        for i in range(0,25):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.MYOFIBROBLAST)
        #for i in range(0,100):
            #randx,randy = randopos(vv.boundaryat, vv.boundaryat, 2)
            #self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.NEUTROPHIL)
        for i in range(0,1000):
            randx,randy = randopos(vv.boundaryat, vv.boundaryat, 3)
            self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.MONOCYTE)
        
        for cell in self.cell_list:
            cell.targetVolume = 1
            cell.lambdaVolume = 10.0
            
            
            if cell.type == 1: # endothelial
                cell.dict["life"] = 0
                cell.dict["span"] = vv.lifespane
                cell.dict["divide"] = vv.divpre
                cell.dict["dividepr"] = 1
                
            if cell.type == 2: # neutrophil
                # print(cell.type)
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL8")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(setSaturationCoef)
                cell.dict["span"] = vv.lifespannr
                #cell.dict["divide"] = random.randint(0,vv.divprnr)
                cell.dict["dividepr"] = vv.divprnr
                cell.dict["life"] = random.randint(0,vv.lifespannr)
                
            if cell.type == 3: # monocyte
                # print(cell.type)
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(setSaturationCoef)
                cell.dict["span"] = vv.lifespanm
                #cell.dict["divide"] = random.randint(0,vv.divprm)
                cell.dict["dividepr"] = vv.divprm
                cell.dict["life"] = random.randint(0,vv.lifespanm)
            if cell.type == 3: # monocyte
                # print(cell.type)
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TNF")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(setSaturationCoef)
                cell.dict["span"] = vv.lifespanm
                #cell.dict["divide"] = random.randint(0,vv.divprm)
                cell.dict["dividepr"] = vv.divprm
                cell.dict["life"] = random.randint(0,vv.lifespanm)
            if cell.type == 4: # fibroblast
                # print(cell.type)
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TGF")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(setSaturationCoef)
                cell.dict["span"] = vv.lifespanf
                cell.dict["divide"] = random.randint(0,vv.divprf)
                cell.dict["dividepr"] = vv.divprf
                cell.dict["life"] = random.randint(0,vv.lifespanf)
            if cell.type == 10:
                # print(cell.type)
                #cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TGF")
                #cd.setLambda(2000.0)
                #cd.setChemotactTowards("MEDIUM")
                #cd.setSaturationCoef(10**-12)
                cell.dict["span"] = vv.lifespanf
                cell.dict["divide"] = random.randint(0,vv.divprf)
                cell.dict["dividepr"] = vv.divprf
                cell.dict["life"] = random.randint(0,vv.lifespanf)
                


    def step(self,mcs):
        
        global cytokines
        global cellpresent
        global fullFileName
        global relaxationmcs

        sigmoida = 1
        sigmoidb = 4
        ccount = np.zeros(vv.total_celltypes+1)
        if mcs%(vv.relaxationmcs*10)==0:
            for i in range(0,10):
                randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
                self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.NEUTROPHIL)
            for i in range(0,10):
                randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
                self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.FIBROBLAST)
            for i in range(0,3):
                randx,randy = randopos(vv.boundaryat, vv.boundaryat, 1)
                self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.MYOFIBROBLAST)
            for i in range(0,100):
                randx,randy = randopos(vv.boundaryat, vv.boundaryat, 3)
                self.cell_field[randx:randx+1, randy:randy+1, 0] = self.new_cell(self.MONOCYTE)
            
            for cell in self.cell_list:
                cell.targetVolume = 1
                cell.lambdaVolume = 10.0
            
                if cell.type == 1: # endothelial
                    cell.dict["life"] = 0
                    cell.dict["span"] = vv.lifespane
                    cell.dict["divide"] = vv.divpre
                    cell.dict["dividepr"] = 1
                    
                if cell.type == 2: # neutrophil
                    # print(cell.type)
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL8")
                    cd.setLambda(setlambda)
                    cd.setChemotactTowards("MEDIUM")
                    cd.setSaturationCoef(setSaturationCoef)
                    cell.dict["span"] = vv.lifespannr
                    #cell.dict["divide"] = random.randint(0,vv.divprnr)
                    cell.dict["dividepr"] = vv.divprnr
                    cell.dict["life"] = random.randint(0,vv.lifespannr)
                    
                if cell.type == 3: # monocyte
                    # print(cell.type)
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1")
                    cd.setLambda(setlambda)
                    cd.setChemotactTowards("MEDIUM")
                    cd.setSaturationCoef(setSaturationCoef)
                    cell.dict["span"] = vv.lifespanm
                    #cell.dict["divide"] = random.randint(0,vv.divprm)
                    cell.dict["dividepr"] = vv.divprm
                    cell.dict["life"] = random.randint(0,vv.lifespanm)
                if cell.type == 3: # monocyte
                    # print(cell.type)
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TNF")
                    cd.setLambda(setlambda)
                    cd.setChemotactTowards("MEDIUM")
                    cd.setSaturationCoef(setSaturationCoef)
                    cell.dict["span"] = vv.lifespanm
                    #cell.dict["divide"] = random.randint(0,vv.divprm)
                    cell.dict["dividepr"] = vv.divprm
                    cell.dict["life"] = random.randint(0,vv.lifespanm)
                if cell.type == 4: # fibroblast
                    # print(cell.type)
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TGF")
                    cd.setLambda(setlambda)
                    cd.setChemotactTowards("MEDIUM")
                    cd.setSaturationCoef(setSaturationCoef)
                    cell.dict["span"] = vv.lifespanf
                    cell.dict["divide"] = random.randint(0,vv.divprf)
                    cell.dict["dividepr"] = vv.divprf
                    cell.dict["life"] = random.randint(0,vv.lifespanf)
                if cell.type == 10:
                    # print(cell.type)
                    #cd = self.chemotaxisPlugin.addChemotaxisData(cell, "TGF")
                    #cd.setLambda(2000.0)
                    #cd.setChemotactTowards("MEDIUM")
                    #cd.setSaturationCoef(10**-12)
                    cell.dict["span"] = vv.lifespanf
                    cell.dict["divide"] = random.randint(0,vv.divprf)
                    cell.dict["dividepr"] = vv.divprf
                    cell.dict["life"] = random.randint(0,vv.lifespanf)

        
        # self.field=CompuCell.getConcentrationField()
        # Make sure PixelTracker plugin is loaded

        # cellsaver = list()
        for cell in self.cell_list:
            xCOM = cell.xCOM
            yCOM = cell.yCOM
            zCOM = cell.zCOM
            if xCOM>=500:
                xCOM = 499
            if yCOM>=500:
                yCOM = 499
            pos=int(xCOM)+(int(yCOM))*nx
            if pos>250000:
                print(xCOM, yCOM)
            #print(cell.type, cell.id)
            cell.dict["life"] += 1
            if cell.type == 1:
                cellpresente[pos]=1.        
            if cell.type == 6:
                cellpresentndn[pos]=1.
            if cell.type == 5:
                cellpresentna[pos]=1.
            if cell.type == 8:
                cellpresentm1[pos]=1.
            if cell.type == 9:
                cellpresentm2[pos]=1.
                
        cytokines = tester(cellpresente,cellpresentndn,cellpresentna,cellpresentm1,cellpresentm2,cytokines, mesh)
        self.scalarFieldil8[:] = np.reshape(cytokines[0], (nx,ny,1), 'F')
        self.scalarFieldil1[:] = np.reshape(cytokines[1], (nx,ny,1), 'F')
        self.scalarFieldil6[:] = np.reshape(cytokines[2], (nx,ny,1), 'F')
        self.scalarFieldil10[:] = np.reshape(cytokines[3], (nx,ny,1), 'F')
        self.scalarFieldtnf[:] = np.reshape(cytokines[4], (nx,ny,1), 'F')
        self.scalarFieldtgf[:] = np.reshape(cytokines[5], (nx,ny,1), 'F')
        
        cellpresente[:] = 0
        cellpresentndn[:] = 0
        cellpresentna[:] = 0
        cellpresentm1[:] = 0
        cellpresentm2[:] = 0
        
        il8_list = []
        il1_list = []
        il6_list = []
        il10_list = []
        tnf_list = []
        tgf_list = []  
        for cell in self.cell_list:
            xCOM = int(cell.xCOM)
            yCOM = int(cell.yCOM)
            zCOM = int(cell.zCOM)
            if xCOM>=500:
                xCOM = 499
            if yCOM>=500:
                yCOM = 499
            ccil8 = self.scalarFieldil8[xCOM,yCOM,zCOM]
            ccil1 = self.scalarFieldil1[xCOM,yCOM,zCOM]
            ccil6 = self.scalarFieldil6[xCOM,yCOM,zCOM]
            ccil10 = self.scalarFieldil10[xCOM,yCOM,zCOM]
            cctnf = self.scalarFieldtnf[xCOM,yCOM,zCOM]
            cctgf = self.scalarFieldtgf[xCOM,yCOM,zCOM]
            
            fileDir=os.path.dirname (os.path.abspath(fullFileName))
            cytoname= fileDir+"/datafiles"+str(mcs)+"concentration.txt"
            if not os.path.exists(cytoname):
                with open(fileDir+"/datafiles"+str(mcs)+"concentration.txt",'w') as cytofile:
                    writer = csv.writer(cytofile)
                    writer.writerow(["mcsteps","xCOM","yCOM","zCOM","il8","il1","il6","il10","tnf","tgf"])
            with open(fileDir+"/datafiles"+str(mcs)+"concentration.txt",'a') as cytofile:
                cytowriter = csv.writer(cytofile)
                cytowriter.writerow([str(mcs),str(xCOM),str(yCOM),str(zCOM),str(ccil8),str(ccil1),str(ccil6),str(ccil10),str(cctnf),str(cctgf)])
                
            if 50<xCOM<450 and 50<yCOM<450:
                il8_list.append(ccil8)
                il1_list.append(ccil1)
                il6_list.append(ccil6)
                il10_list.append(ccil10)
                tnf_list.append(cctnf)
                tgf_list.append(cctgf)
            
            if cell.dict["life"] > cell.dict["span"]:
                self.delete_cell(cell)
                
            else:
                if cell.dict["life"]>vv.timeforgrowth*cell.dict["span"]:
                    cell.targetVolume = 2
                if cell.type == 2: # 2 neutrophil-->5 neutrophila
                    # proba = ccil8/(ccil8 + vv.cil8)*vv.lnril8 + ccil6/(ccil6 + vv.cil6)*vv.lnril6+ ccil1/(ccil1 + vv.cil1)*vv.lnril1+ cctnf/(cctnf + vv.ctnf)*vv.lnrtnf+ccil10/(ccil10 + vv.cil10)*vv.tnril10
                    proba = sigmoid(ccil8*10**9, sigmoida, sigmoidb)*vv.lnril8 + sigmoid(ccil6*10**11, sigmoida, sigmoidb)*vv.lnril6 + sigmoid(ccil1*10**9, sigmoida, sigmoidb)*vv.lnril1 + sigmoid(cctnf*10**9, sigmoida, sigmoidb)*vv.lnrtnf - sigmoid(ccil10*10**12, sigmoida, sigmoidb)*vv.tnril10
                    # print (proba, cell.type)
                    if proba > random.random():
                        cell.type = 5
                if cell.type == 7:  # 7 monocyter-->8 macrophage1
                    # proba = ccil6/(ccil6 + vv.cil6)*vv.lmril6 + cctnf/(cctnf + vv.ctnf)*vv.lmrtnf + ccil10/(ccil10 + vv.cil10)*vv.tmril10
                    proba = 0.1 + 0.9*sigmoid(ccil6*10**11, sigmoida, sigmoidb)*vv.lmril6 + sigmoid(cctnf*10**9, sigmoida, sigmoidb)*vv.lmrtnf - sigmoid(ccil10*10**12, sigmoida, sigmoidb)*vv.tmril10
                    # print (proba, cell.type)
                    if proba > random.random():
                        cell.type = 8
                        cell.dict["life"] = 0 
                        cell.dict["span"] = vv.lifespanmr
                if cell.type == 8:  # 8 macrophage1--> 9 macrophage2 (have base prob)
                    # proba = ccil10/(ccil10 + vv.cil10)*vv.lm1il10
                    # print (proba, cell.type)
                    proba = 0.1 + 0.9*sigmoid(ccil10*10**12, sigmoida, sigmoidb)*vv.lm1il10
                    if proba > random.random():
                        cell.type = 9
                        cell.dict["life"] = 0 
                        cell.dict["span"] = modmm(mcs, vv.km, vv.vmax) #vv.lifespanmr         
                if cell.type == 4:#4 fibroblast-->10 myofibroblast
                    # proba = cctgf/(cctgf + vv.ctgf)*vv.lftgf
                    # print (proba, cell.type)
                    proba = 0.1*sigmoid(cctgf*10**10, sigmoida, sigmoidb)*vv.lftgf
                    if proba > random.random():
                        cell.type = 10
                if cell.type == 5: # have to fix properly 5 neutrophlia-->6 neutrophilndn
                    # print(cell.dict["divide"])
                    #if random.randint(0,int(cell.dict["dividepr"])) == cell.dict["divide"]:
                    # 100/1000 chance necrosis could change 
                    if random.randint(0,1000) == 50:
                        cell.type = 6
                if cell.type == 3: # 3 monocyte--> 7 monocyter (have base prob)
                    # print(ccil8, vv.cil8) because of il6
                    proba = 0.1 + 0.9*sigmoid(ccil6*10**11, sigmoida, sigmoidb)*vv.tranril6
                    if proba > random.random(): 
                        cell.type = 7
                        cell.dict["life"] = 0 
                        cell.dict["span"] = vv.lifespanmr
                ccount[cell.type]+=1
        # print (ccount)
        # save ccount
        self.plot_win.add_data_point('Endothelial', mcs, ccount[1])
        self.plot_win.add_data_point('Neutrophils', mcs, ccount[2])
        self.plot_win.add_data_point('Monocytes', mcs, ccount[3])
        self.plot_win.add_data_point('Fibroblast', mcs, ccount[4])
        self.plot_win.add_data_point('Neutrophil A', mcs, ccount[5])
        self.plot_win.add_data_point('Neutrophil NDN', mcs, ccount[6])
        self.plot_win.add_data_point('Monocyte R', mcs, ccount[7])
        self.plot_win.add_data_point('Macrophage I', mcs, ccount[8])
        self.plot_win.add_data_point('Macrophage II', mcs, ccount[9])
        self.plot_win.add_data_point('Myofibroblast', mcs, ccount[10])
        fileDir=os.path.dirname (os.path.abspath(fullFileName))
        namer=fileDir+"/datafiles"+str(mcs)+".png"
        if not os.path.exists(os.path.dirname(namer)):
            os.makedirs(os.path.dirname(namer))
        self.plot_win.save_plot_as_png(namer, 1200, 1200)

        countname = fileDir+"/cellcount.txt"
        if not os.path.exists(countname):
            with open(fileDir+"/cellcount.txt",'w') as file:
                writer = csv.writer(file)
                writer.writerow(["mcsteps","1","2","3","4","5","6","7","8","9","10"])
        with open(fileDir+"/cellcount.txt", "a") as file:
            writer = csv.writer(file)
            writer.writerow([str(mcs),str(ccount[1]),str(ccount[2]),str(ccount[3]),str(ccount[4]),str(ccount[5]),str(ccount[6]),str(ccount[7]),str(ccount[8]),str(ccount[9]),str(ccount[10])])
        
        il8_mean = np.mean(il8_list)
        il1_mean = np.mean(il1_list)
        il6_mean = np.mean(il6_list)
        il10_mean = np.mean(il10_list)
        tnf_mean = np.mean(tnf_list)
        tgf_mean = np.mean(tgf_list)
        
        il8_std = np.std(il8_list)
        il1_std = np.std(il1_list)
        il6_std = np.std(il6_list)
        il10_std = np.std(il10_list)
        tnf_std = np.std(tnf_list)
        tgf_std = np.std(tgf_list)
            
        cytodata = fileDir+"/mean_concentration.txt"
        if not os.path.exists(cytodata):
            with open(fileDir+"/mean_concentration.txt",'w') as cytofile:
                cytowriter = csv.writer(cytofile)
                cytowriter.writerow(["meanconcen","il8mean","il1mean","il6mean","il10mean","tnfmean","tgfmean","il8std","il1std","il6std","il10std","tnfstd","tgfstd"])
        with open(fileDir+"/mean_concentration.txt", "a") as cytofile:
            cytowriter = csv.writer(cytofile)
            cytowriter.writerow([str(mcs),str(il8_mean),str(il1_mean),str(il6_mean),str(il10_mean),str(tnf_mean),str(tgf_mean),str(il8_std),str(il1_std),str(il6_std),str(il10_std),str(tnf_std),str(tgf_std)])
        
        
        
            
            
    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return


        
class celldivisionSteppable(MitosisSteppableBase):
    def __init__(self,frequency=int(vv.relaxationmcs)):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list:
            # about cell division!!
            # if cell.volume>=2 and cell.dict["life"] > vv.timeforgrowth*cell.dict["span"] and random.randint(0,int(cell.dict["dividepr"])) == cell.dict["divide"]:
            if cell.volume>=2 and cell.dict["life"] > vv.timeforgrowth*cell.dict["span"] and random.randint(0,1000) == 5000:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:
            self.divide_cell_along_major_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume = 1.0                
        self.parent_cell.dict["life"] = 0
        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        

        
