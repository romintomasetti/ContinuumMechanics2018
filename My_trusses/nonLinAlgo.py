import sys, os
import numpy as np
from truss import *
import matplotlib.pyplot as plt

class NonLinearAlgorithm: # Father class for all the non-linear algorithms
    
    def __init__(self, _truss, _toll, _nItMax,_computeTangentMethod):
        self.truss = _truss                     # The truss
        self.toll = _toll                       # Tolerance used to assess convergence
        self.nItMax = _nItMax                   # Maximum number of iterations
        self.lamdaMax = 1.0                     # Maximum value for lamda (a priori, it can be different from 1)
        self.stopit = False                     # Only used by the GUI
        self.cleanWorkspace()                   # Run at the beginning of the simulation to clean the workspace: all the old result files will be deleted
        self.computeTangentMethod = _computeTangentMethod  # Or 'numeric' -> Used in bar.buildKt()
        self.dhook = None                       # Only used by the GUI
        
    def display(self, step, lamda):
        if self.dhook:
            self.dhook.display(step, lamda)
            self.dhook.refresh()
        
    def archive(self, lamda, nIt): # Writes the files containing the results. Just the main information is archived, but you can enrich it if you want...
        f1 = open(('Lambda.ascii'),'a')
        f1.write(str(lamda)+'\n')
        f1.close()
        f2 = open(('Iterations.ascii'),'a')
        f2.write(str(nIt)+'\n')
        f2.close()
        for node in self.truss.nodes:
            f3 = open(('Node_'+str(node.nb)+'_POSITION.ascii'),'a')
            f3.write(str(node.x)+'   '+str(node.y)+'\n')
            f3.close()
            f4 = open(('Node_'+str(node.nb)+'_DISPLACEMENTS.ascii'),'a')
            f4.write(str(node.u)+'   '+str(node.v)+'\n')
            f4.close()
            f5 = open(('Node_'+str(node.nb)+'_F_EXT.ascii'),'a')
            f5.write(str(node.Fx)+'   '+str(node.Fy)+'\n')
            f5.close()
            
    def archive_internal_forces(self):
        #Stores the internal forces inside each bar:
        for bar in self.truss.bars:
            f1 = open(('Bar_'+str(bar.nb)+'_internal_force.ascii'),'a')
            i_f = bar.getFint()
            f1.write(str(i_f[0])+'   '+str(i_f[1])+'   '+str(i_f[2])+'   '+str(i_f[3])+'\n')
            f1.close()
    
    def computeError(self,g,lambda_): # Computes the error used to assess convergence
        error = np.linalg.norm(self.truss.getOOBF())*self.toll
        return error
    
    def cleanWorkspace(self):
        dir_name = os.getcwd()
        test = os.listdir(dir_name)
        
        for item in test:
            if item.endswith(".ascii"):
                os.remove(os.path.join(dir_name, item))
    
class IncrementalAlgorithm(NonLinearAlgorithm): # A simple explicit solver for non-linear problems
    
    def __init__(self, _truss, _dlamda): 
        
        NonLinearAlgorithm.__init__(self, _truss, 0., 0)
        self.dlamda = _dlamda
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0
        while lamda < self.lamdaMax and not self.stopit:
            step+=1
            lamda = lamda0 + dlamda
            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            
            print '\n--- Incremental method - Load level, lambda =', lamda, ' ---'
            
            self.truss.applyNodalLoads(lamda)
            
            Kt = self.truss.buildKt(self.computeTangentMethod)
            g = self.truss.getOOBF()
            du = np.linalg.solve(Kt, -g)
            
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            
            self.truss.update()
            self.archive(lamda, 0)
            self.display(step, lamda)
            lamda0 = lamda

class NewtonRaphsonAlgorithm(NonLinearAlgorithm): # Newton-Raphson method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax)
        self.dlamda = _dlamda
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0
        print ""
        while lamda < self.lamdaMax and not self.stopit:
            step+=1
            lamda = lamda0 + dlamda
            Delta_u = np.zeros((2*len(self.truss.nodes),1)) 
            
            print '>>> Newton Raphson - Load level, lambda = ', lamda, ' --- dlambda = ',dlamda,'        \r',
            
            #Apply loads:
            self.truss.applyNodalLoads(lamda)
            
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Calcule la difference entre les forces internes et externes
            g = self.truss.getOOBF()  
            # solve a linear matrix equation
            du = np.linalg.solve(Kt, -g)  
                
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            #Compute the error:
            error = self.computeError(g, lamda)
            #Apply the corrector phase until tolerance is satisfied:
            current_iteration = 0
            while error > self.toll and current_iteration < self.nItMax:
                #Build the stiffness matrix                
                Kt = self.truss.buildKt(self.computeTangentMethod)
                #Compute the out of balance forces:
                g_oof = self.truss.getOOBF()
                #Solve the system:
                du = np.linalg.solve(Kt, -g_oof)
                
                Delta_u+=du
                self.truss.incrementPositions(Delta_u)
                error = self.computeError(g, lamda)
                current_iteration +=1
                
            self.truss.update()
            self.archive(lamda, current_iteration)
            self.display(step, lamda)
            lamda0 = lamda
        return 1

class QuasiNewtonAlgorithm(NonLinearAlgorithm): # Quasi-Newton method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _Kt_method): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax,_Kt_method)
        self.dlamda = _dlamda
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0
        while lamda < self.lamdaMax and not self.stopit:
            step+=1
            lamda = lamda0 + dlamda
            Delta_u = np.zeros((2*len(self.truss.nodes),1)) #*2 pour les 2 directions x et y 
            
            print '\n--- Newton Raphson - Load level, lambda =', lamda, ' ---'
            
            self.truss.applyNodalLoads(lamda)
            
            Kt = self.truss.buildKt(self.computeTangentMethod)
            g  = self.truss.getOOBF()
            try:
                Kt_inv = np.linalg.inv(Kt)
            except np.linalg.LinAlgError:
                print "Matrix Kt not invertible"
                print(Kt)
                sys.exit()
            else:
                du = np.dot(Kt_inv, -g)
                Delta_u+=du
                self.truss.incrementPositions(Delta_u)
                error = 1.0
                error = self.computeError(g, lamda)
                NbIter = 0
                while error > self.toll and NbIter < self.nItMax:                
                    g = self.truss.getOOBF()
                    du = np.dot(Kt_inv, -g) 
                
                    Delta_u+=du
                    self.truss.incrementPositions(Delta_u)
                    error = self.computeError(g, lamda)
                    NbIter +=1
                
                self.truss.update()
                self.archive(lamda, NbIter)
                self.display(step, lamda)
                lamda0 = lamda
        return 1
        
class NewtonRaphsonAlgorithm_inTime(NonLinearAlgorithm):
    #Constructor:
    def __init__(self, _truss, _toll, _nItMax, _dt, _t_max,_F_max, _period):
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax,'')
        #Time step:
        self.dt = _dt
        #Maximum time before stopping:
        self.t_max = _t_max
        #Maximum force:
        self.F_max = _F_max
        #Period of the signal:
        self.period = _period
        
    #Main function, only valid for a truss with 2 nodes !
    #All loads are applied on the node 2 ! In the Y direction.
    # DDL(X,Y) are fixed for node 1.
    # DDL(X)   is  fixed for node 2.
    def run(self):
        #Check the time step that is given to be sure we go through all minimas and maximas
        #    of the sawtooth loading:
        if not self.period/self.dt == math.floor(self.period/self.dt):
            sys.exit('The provided time step is not correct.')
        
        #Initialize the time:
        t = 0.
        
        #Initialize nodal loads:
        self.truss.nodes[1].applyTimeLoad('y', self.F_sawtooth)
        
        dt = self.dt
        step = 0
        #Loop on time:
        while t < self.t_max:
            #Increment the step:
            step+=1
            
            #Increment time:
            t+=dt
            
            #Allocated space:
            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            
            print '>>> Newton Raphson in time - Load time, t =', t, ' ---         \r',
            
            #Apply loads:
            self.truss.applyNodalTimeLoads(t)
            #Build the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Compute the out of balance forces:
            g = self.truss.getOOBF()
            #Solve the linear system:
            du = np.linalg.solve(Kt, -g)
            #Update positions:
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            #Compute the error:
            error = self.computeError_time(g, t)
            
            #Usual N.-R. step:
            current_iteration = 0
            while error > self.toll :
                if current_iteration > self.nItMax:    
                    print ""
                    print ">>> N.-R. in time failed at converging."
                    sys.exit()
                #Compute the tangent stiffness matrix:
                Kt = self.truss.buildKt(self.computeTangentMethod)
                #Compute the out of balance forces:
                g = self.truss.getOOBF()
                #Solve the linear system:
                du = np.linalg.solve(Kt, -g)
                #Update positions:
                Delta_u+=du
                self.truss.incrementPositions(Delta_u)
                #Compute the error:
                error = self.computeError_time(g, t)
                current_iteration += 1
            
            self.truss.update()
            self.archive(t, current_iteration)
            self.archive_some_time_infos(t)
            self.archive_PK2(t)
            self.display(step, t)
        
        return 1
        
    #Archive Pk2 stress in the bars:
    def archive_PK2(self,t):
        f1 = open('NR_in_time_PK2.ascii','a')
        f1.write('t=' + str(t))
        for bar in self.truss.bars:
            pk2 = bar.getPK2Stress()
            f1.write(';PK2_bar' + str(bar.nb) + '=' + str(pk2))
        f1.write('\n')
        f1.close()
        
    #Compute a saw-tooth-like loading:
    def F_sawtooth(self, t):
        t_equ = t/self.period - math.floor(t/self.period) # t_equ is the equivalent time which has the same force value.
        if(t_equ <= 0.25):
            F = 4.0*self.F_max*t_equ
        elif (t_equ > 0.25 and t_equ <= 0.75):
            F = 2.0*self.F_max*(-2.0*t_equ+1.0)
        else:
            F = 4.0*self.F_max*(t_equ-1.0)
        return F
        
    def archive_some_time_infos(self,t):
        f1 = open('Some_more_infos.ascii','a')
        f1.write('t=' + str(t) + ';dt=' + str(self.dt) + ';F_max=' + str(self.F_max))
        f1.write(';T=' + str(self.period) + '\r')
        f1.close()
    
    def computeError_time(self, g, t):
        
        error = abs(np.linalg.norm(g))*self.toll

        return error

class ArcLengthAlgorithm(NonLinearAlgorithm): # Arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id, _applyCorrectiveMethod,_computeTangentMethod): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax,_computeTangentMethod)
        self.dlamda = _dlamda
        self.psi = _psi
        self.Id = _Id # Ideal number of iterations per step
        self.nItMax = _nItMax
        self.applyCorrectiveMethod = _applyCorrectiveMethod
        
    def archive_truss_(self,step,lamda):
            str_ = '_'+str(self.dlamda)+'_'+str(self.Id)+'_'+str(self.psi)
            str_ = str_+'_'+'_'+str(self.nItMax)+'_'+str(self.toll)
            
            f1 = open(('Lambda'+str_+'.ascii'),'a')
            f1.write(str(lamda)+'\n')
            f1.close()
            for node in self.truss.nodes:
                f3 = open(('Node_'+str(node.nb)+'_POSITION'+str_+'.ascii'),'a')
                f3.write(str(node.x)+'   '+str(node.y)+'\n')
                f3.close()
                f4 = open(('Node_'+str(node.nb)+'_DISPLACEMENTS'+str_+'.ascii'),'a')
                f4.write(str(node.u)+'   '+str(node.v)+'\n')
                f4.close()
                f5 = open(('Node_'+str(node.nb)+'_F_EXT'+str_+'.ascii'),'a')
                f5.write(str(node.Fx)+'   '+str(node.Fy)+'\n')
                f5.close()
        
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0    
        DeltaL = 0.
        DeltaL0 = 0.
        current_iteration = 0
        # The first step is a call to the Newton-Raphson algorithm:
        
        step+=1
        lamda = lamda0 + dlamda
        Delta_u = np.zeros((2*len(self.truss.nodes),1)) 
        
        print '\n--- ArcLengthAlgorithm :: Newton-Raphson step - Load level, lambda =', lamda, ' ---'
        
        #Apply loads:
        self.truss.applyNodalLoads(lamda)
        #Build the tangent stiffness matrix:
        Kt = self.truss.buildKt(self.computeTangentMethod)
        #Get the out of balance forces:
        g = self.truss.getOOBF()
        #Solve the linear system for du:
        du = np.linalg.solve(Kt, -g)
        #Assemble the global loads vector:
        qef = self.truss.get_qef()
        #Increment du:
        Delta_u+=du
        #Update positions:
        self.truss.incrementPositions(Delta_u)
        #Compute the error:
        error = self.computeError(g, lamda)
        #While loop of the Newton-Raphson step:
        while error > self.toll and current_iteration < self.nItMax:         
            Kt = self.truss.buildKt(self.computeTangentMethod)
            g = self.truss.getOOBF()
            du = np.linalg.solve(Kt, -g)            
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            error = self.computeError(g, lamda)
            current_iteration +=1
        if current_iteration > self.nItMax:
            print ">>> N.-R. step didn't converge in the spherical arc-length method."
            sys.exit()
        #Get the initial arc-length:
        DeltaL0 = math.sqrt(np.transpose(Delta_u).dot(Delta_u) + self.psi**2*lamda**2*np.transpose(qef).dot(qef)) # !!!!!!!!!!!!!! A refaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DeltaL = DeltaL0
        lamda0 = lamda
        #Update the truss:
        self.truss.update()
        self.archive(lamda, current_iteration)
        self.display(step, lamda) 
        self.archive_truss_(step,lamda)
        
        #########################################
        ## PREDICTOR PHASE ##
        #####################
        
        #We want to track oscillations:
        counter_oscillating = 0
        boolean_oscillating = 0
        counter_applied_corrective = 0
        val = 0
        bool_corr = 0
        
        RESTART_ = 0
        
        while lamda < self.lamdaMax and not self.stopit:
            if step == 6e3:
                return -7
            
            #Not yet implemented !
            if RESTART_ == 1:
                self.truss.resetPositions()
                print ''
                print 'truss has been reset to previous !',DeltaL,current_iteration
            
            dlamda_previous = dlamda
            
            step+=1
            #Apply loads:
            self.truss.applyNodalLoads(lamda)
     
            current_iteration = 0

            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            #Build the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Get the out of balance forces:
            g = self.truss.getOOBF()
            #Get the global loads vector:
            qef = self.truss.get_qef() 
            #Try to invert Kt:
            try:
                Kt_inv = np.linalg.inv(Kt)
            except np.linalg.LinAlgError:
                # Not invertible. Skip this one.
                print "Kt is not invertible"
                print(Kt)
                sys.exit()
            else:
                deltaP = - Kt_inv.dot(g) 
                deltaPt = Kt_inv.dot(qef)
                #Compute deltaLambda_p:
                Ct1 = np.transpose(deltaP).dot(deltaP)*np.transpose(deltaPt).dot(deltaPt)+ (DeltaL**2-np.transpose(deltaP).dot(deltaP))*(np.transpose(deltaPt).dot(deltaPt)+self.psi**2*np.transpose(qef).dot(qef))
                if Ct1 < 0:
                    print "Computation of Ct1 leads to a negative number. We can't compute its square root."
                    sys.exit()
                else:
                    Ct1 = math.sqrt(Ct1)

                # Verify positive-definiteness of the tangent stiffness matrix:
                if np.all(np.linalg.eigvals(Kt) > 0):
                    #Take the positive sign '+':
                    dlamda = (-np.transpose(deltaP).dot(deltaPt) + Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))       
                else:
                    dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))       
                   
                ##############################
                ## APPLY CORRECTIVE PROCESS ##
                ##############################
                MAXX = 300
                if self.applyCorrectiveMethod == 1 and counter_oscillating > 10 and counter_applied_corrective < MAXX:
                    counter_applied_corrective += 1
                    if counter_applied_corrective == 1:
                        du_corr = np.linalg.solve(Kt_inv, -g);
                        val = qef.T.dot(du_corr)
                    if val < 0:
                        dlamda = (-np.transpose(deltaP).dot(deltaPt) + Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                    else:
                        dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                    bool_corr = 1
                
                #Update du:
                du = deltaP + dlamda*deltaPt
                Delta_u+=du
                lamda = lamda0 + dlamda
                self.truss.incrementPositions(Delta_u)
                error = self.computeError(g, Kt)
                
                #####################
                ## CORRECTOR PHASE ##
                #####################
                
                while error > self.toll and current_iteration < self.nItMax:
                    #Apply loads:
                    self.truss.applyNodalLoads(lamda)
                    #Build the tangent stiffness matrix:
                    Kt = self.truss.buildKt(self.computeTangentMethod)
                    #Get the initial out of balance forces:
                    g0 = self.truss.getOOBF()
                    #Get the globl loads vector:
                    qef = self.truss.get_qef()
                    try:
                        Kt_inv = np.linalg.inv(Kt)
                    except numpy.linalg.LinAlgError:
                        print "Kt is not invertible"
                        printf(Kt)                        
                        sys.exit()
                    else:
                        deltaP = - Kt_inv.dot(g0) 
                        deltaPt =  Kt_inv.dot(qef)
                        #Compute coefficients of the equation (15):
                        a1 = np.transpose(deltaPt).dot(deltaPt) + self.psi**2*np.transpose(qef).dot(qef)
                        a2 = 2 * (np.transpose(deltaPt).dot(Delta_u+deltaP) + dlamda*self.psi**2*np.transpose(qef).dot(qef))
                        a3 = np.transpose(Delta_u+deltaP).dot(Delta_u+deltaP) + dlamda**2 *self.psi**2*np.transpose(qef).dot(qef) - DeltaL**2
                        #Compute the roots of the second order equation (15):
                        try:
                            lamda1 = ( -a2 + math.sqrt(a2**2 - 4*a1*a3)) / (2*a1)
                            lamda2 = ( -a2 - math.sqrt(a2**2 - 4*a1*a3)) / (2*a1)
                        except:
                            print a1 , a2 , a3, a2**2 - 4*a1*a3, lamda
                            sys.exit()
                        du1 = deltaP + lamda1*deltaPt
                        du2 = deltaP + lamda2*deltaPt
                        #Compute the cosinus of the angle:
                        cTheta1 = (np.transpose(Delta_u).dot(Delta_u+du1) + self.psi**2 * np.transpose(qef).dot(qef) * dlamda * (lamda1+dlamda))                        
                        cTheta2 = (np.transpose(Delta_u).dot(Delta_u+du2) + self.psi**2 * np.transpose(qef).dot(qef) * dlamda * (lamda2+dlamda))
                        #Choose the smallest:
                        if cTheta1>cTheta2:
                            dlamda += lamda1
                            du = du1
                        else:
                            dlamda += lamda2
                            du = du2
                        #Update:
                        Delta_u += du
                        lamda = lamda0 + dlamda 
                        self.truss.incrementPositions(Delta_u)
                        #Compute the error:
                        error = self.computeError(g0, lamda)
                        current_iteration += 1
                        
                # if error > self.toll:
                    # RESTART_ = 1
                # else:
                    # RESTART_ = 0
                       
                if current_iteration != 0 :#and RESTART_ == 0:    
                    #Update deltaL:
                    DeltaL = DeltaL0*math.sqrt(float(self.Id)/current_iteration) 

                else:#if RESTART_ == 0:
                    #update the deltaL variable with a default step:
                    DeltaL = DeltaL0*math.sqrt(float(self.Id)/1.2)
                    #Note: 1.2 is chosen quite arbitrarily.
                # else:
                    # DeltaL = 0.8*DeltaL
                
                if RESTART_ == 0:
                    print ">> Sph. arc-len. :: pred. at step ",step," ----- lambda: ",lamda," --- dlambda: ",dlamda,"--- dlambda_previous: ",dlamda_previous,"         \r",
                    #Update the truss:
                    self.truss.update()
                    self.archive(lamda[0,0], current_iteration)
                    self.display(step, lamda)
                    self.archive_truss_(step,lamda)
                    #self.archive_internal_forces()
                    lamda0 = lamda
                
                #Check that we don't oscillate between two values of dlamda:
                if abs(abs(dlamda_previous)-abs(dlamda)) < 1e-5:
                    #print "Oscillating!"
                    counter_oscillating += 1
                    # Check the signs were opposite!
                    if dlamda_previous < 0 and dlamda < 0:
                        counter_oscillating -= 1
                    if dlamda_previous > 0 and dlamda > 0:
                        counter_oscillating -= 1
                    if counter_oscillating > 50 and not self.applyCorrectiveMethod == 1:
                        #boolean_oscillating = 1
                        counter_oscillating = 0
                        print ""
                        print ">>> Oscillations."
                        return -77
                else:
                    #print " >>> In undef.  : ",abs(abs(dlamda_previous_at_beginning_of_pred)-abs(dlamda))
                    if not self.applyCorrectiveMethod == 1:
                        counter_oscillating = 0
                        boolean_oscillating = 0
                    
                if step > 10000 and self.applyCorrectiveMethod == 1:
                    sys.exit("Je quitte a 10000 !")
        
        return 1
        

class UpdatedNormalPlaneArcLengthAlgorithm(ArcLengthAlgorithm): # Updated normal plane arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id):
        ArcLengthAlgorithm.__init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id)
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0    
        DeltaL = 0.
        DeltaL0 = 0.
        current_iteration = 0
        #################################
        ## FIRST STEP: N.-R. ALGORITHM ##
        #################################
        
        step+=1
        lamda = lamda0 + dlamda
        Delta_u = np.zeros((2*len(self.truss.nodes),1))
        
        #print '\n---UpdatedNormalPlane :: Newton Raphson - Load level, lambda =', lamda, ' ---'
        
        #Apply loads:
        self.truss.applyNodalLoads(lamda)
        #Build the tangent stiffness matrix:
        Kt = self.truss.buildKt(self.computeTangentMethod)
        #Compute the out of balance forces:
        g = self.truss.getOOBF()
        #Solve the linear system:
        du = np.linalg.solve(Kt, -g)
        #Get the global loads vector:
        qef = self.truss.get_qef()
            
        Delta_u+=du
        #Update positions:
        self.truss.incrementPositions(Delta_u)
        #Compute the error:
        error = self.computeError(g, lamda)
        current_iteration = 0
        #Newton-Raphson loop:
        while error > self.toll and current_iteration < self.nItMax:
            #Compute the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Get the out of balance forces:
            g = self.truss.getOOBF()
            #Solve the system:
            du = np.linalg.solve(Kt, -g)
            #Update positions:
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            #Compute the error:
            error = self.computeError(g, lamda)
            current_iteration +=1
        if current_iteration > self.nItMax:
            print "Update normal plane::Newton-Raphson step didn't converge."
            return -7
            #sys.exit()
        #Compute deltaL:
        DeltaL0 = math.sqrt(np.transpose(Delta_u).dot(Delta_u) + self.psi**2*lamda**2*np.transpose(qef).dot(qef)) # !!!!!!!!!!!!!! A refaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DeltaL = DeltaL0
        lamda0 = lamda
        #Update the truss:
        self.truss.update()
        self.archive(lamda,current_iteration)
        self.display(step, lamda) 
        
        boolean_oscillating = 0
        counter_oscillating = 0
        
        ####################
        ## PREDICTOR STEP ##
        ####################
        while lamda < self.lamdaMax and not self.stopit:
            dlamda_previous_at_beginning_of_pred = dlamda
                        
            step+=1
            #Apply loads:
            self.truss.applyNodalLoads(lamda)
            
            current_iteration = 0
            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            #Compute the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Get the out of balance forces:
            g = self.truss.getOOBF()
            qef = self.truss.get_qef()
            #Try to invert the tangent stiffness matrix:
            try:
                Kt_inv = np.linalg.inv(Kt)
            except numpy.linalg.LinAlgError:
                # Not invertible. Skip this one.
                #print "Kt is not invertible"
                print(Kt)
                #sys.exit()
                return -7
            else:
                deltaP = - Kt_inv.dot(g)
                deltaPt = Kt_inv.dot(qef)
                
                Ct1 = np.transpose(deltaP).dot(deltaP)*np.transpose(deltaPt).dot(deltaPt)+ (DeltaL**2-np.transpose(deltaP).dot(deltaP))*(np.transpose(deltaPt).dot(deltaPt)+self.psi**2*np.transpose(qef).dot(qef))
                if Ct1 < 0:
                    print "Ct1 is negative. Can't take its square root."
                    return -7
                    #sys.exit()
                else:
                    Ct1 = math.sqrt(Ct1)
                try:
                    if np.all(np.linalg.eigvals(Kt) > 0) or (step > 45 and step < 60) :
                        if boolean_oscillating == 0:
                            dlamda = (-np.transpose(deltaP).dot(deltaPt) + Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                            #print "def. pos. -- usual"
                        else:
                            dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                            #print "def. pos. -- changed"
                    elif np.all(np.linalg.eigvals(Kt) < 0):
                        dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                        print "Def. neg."
                        return -7
                        #sys.exit()
                    else:
                        if boolean_oscillating == 0:
                            dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                            #print "undef. -- usual"
                        else:
                            dlamda = (-np.transpose(deltaP).dot(deltaPt) + Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                            #print "undef. -- changed"
                except:
                    print(Kt)
                    fig = plt.figure()
                    plt.imshow(Kt,interpolation='none')
                    plt.show()
                    print "Unexpected error:", sys.exc_info()[0]
                    #sys.exit()
                    return -7
                
                du = deltaP + dlamda*deltaPt
                Delta_u+=du
                lamda = lamda0 + dlamda
                #Increment position:
                self.truss.incrementPositions(Delta_u)
                #Compute error:
                error = self.computeError(g, Kt)
                
                #####################
                ## CORRECTOR PHASE ##
                #####################
                while error > self.toll and current_iteration < self.nItMax:
                    #Apply loads:
                    #print "Corrector",du
                    self.truss.applyNodalLoads(lamda)
                    #Compute the tangent stiffness matrix:
                    Kt = self.truss.buildKt(self.computeTangentMethod)
                    g0 = self.truss.getOOBF()
                    qef = self.truss.get_qef()
                    try:
                        Kt_inv = np.linalg.inv(Kt)
                    except numpy.linalg.LinAlgError:
                        print "Kt is not invertible"
                        print(Kt)
                        print "Unexpected error:", sys.exc_info()[0]
                        #sys.exit()
                        return -7
                    else:
                        deltaP  = - Kt_inv.dot(g0)
                        deltaPt =  Kt_inv.dot(qef)

                        a1 = dlamda*self.psi**2*np.transpose(qef).dot(qef)
                        a1 = a1 + np.transpose(deltaPt).dot(Delta_u)
                        a2 = np.transpose(Delta_u).dot(deltaP)
                        
                        #print "a1 ",a1,"a2 ",a2,"ratio ",a2/a1,"deltzPt ",deltaPt

                        lamda1 = -a2/a1
                        dlamda += lamda1

                        du = deltaP + lamda1*deltaPt

                        Delta_u += du
                        lamda = lamda0 + dlamda
                        #Increment positions:
                        self.truss.incrementPositions(Delta_u)
                        #Compute the error:
                        error = self.computeError(g0, lamda)
                        current_iteration += 1
                        
                       
                if current_iteration != 0:     
                    DeltaL = DeltaL0*math.sqrt(float(self.Id)/current_iteration)

                else :
                    print "Don't go over here !!"
                    sys.exit("Don't go over here!")
                    return -7
                    #sys.exit()
                        
                    
                self.truss.update()
                self.archive(lamda[0,0], current_iteration)
                self.display(step, lamda)
                lamda0 = lamda
                print "- Update normal plane :: The predictor at step ",step," ------- lambda: ",lamda,"---- dlambda: ",dlamda," -- dlamda_prev: ",dlamda_previous_at_beginning_of_pred,"            \r",
                
                #Check that we don't oscillate between two values of dlamda:
                if abs(abs(dlamda_previous_at_beginning_of_pred)-abs(dlamda)) < 1e-6:
                    #print "Oscillating!"
                    counter_oscillating += 1
                    # Check the signs were opposite!
                    if dlamda_previous_at_beginning_of_pred < 0 and dlamda < 0:
                        counter_oscillating -= 1
                    if dlamda_previous_at_beginning_of_pred > 0 and dlamda > 0:
                        counter_oscillating -= 1
                    if counter_oscillating > 50:
                        #boolean_oscillating = 1
                        counter_oscillating = 0
                        print ""
                        print ">>> Oscillations."
                        return -7
                else:
                    #print " >>> In undef.  : ",abs(abs(dlamda_previous_at_beginning_of_pred)-abs(dlamda))
                    counter_oscillating = 0
                    boolean_oscillating = 0

        print "Lambda avant le return: ",lamda                
        return 1

class NormalPlaneArcLengthAlgorithm(ArcLengthAlgorithm): # Normal plane arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id):
        ArcLengthAlgorithm.__init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id)
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0    
        DeltaL = 0.
        DeltaL0 = 0.
        current_iteration = 0
        dLamdaPredict = 0.0
        #################################
        ## FIRST STEP : NEWTON-RAPHSON ##
        #################################
        step+=1
        lamda = lamda0 + dlamda
        Delta_u = np.zeros((2*len(self.truss.nodes),1)) 
        
        print '\n---NormalPlaneArcLengthAlgorithm :: Newton Raphson - Load level, lambda =', lamda, ' ---'
        
        #Apply loads:
        self.truss.applyNodalLoads(lamda)
        #Compute tangent stiffness matrix:
        Kt = self.truss.buildKt(self.computeTangentMethod)
        #Get the out of balance forces:
        g = self.truss.getOOBF()  
        #Solve the linear system:
        du = np.linalg.solve(Kt, -g)  
        #Get the global loads vector:
        qef = self.truss.get_qef()
        #Update deltaU:
        Delta_u+=du
        #Update positions:
        self.truss.incrementPositions(Delta_u)
        #Compute the error:
        error = self.computeError(g, lamda)
        current_iteration = 0
        #Loop of the Newton-Raphson algorithm:
        while error > self.toll and current_iteration < self.nItMax:
            #Get the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Get the out of balance forces:
            g = self.truss.getOOBF()
            #Solve the linear system:
            du = np.linalg.solve(Kt, -g) 
            #Update du:
            Delta_u+=du
            #Update positions:
            self.truss.incrementPositions(Delta_u)
            #Compute the error:
            error = self.computeError(g, lamda)
            current_iteration +=1
        if current_iteration > self.nItMax:
            print "NormalPlaneArcLengthAlgorithm :: N.-R. didn't converge."
            sys.out()
        #Compute deltaL_0 and deltaL:
        DeltaL0 = math.sqrt(np.transpose(Delta_u).dot(Delta_u) + self.psi**2*lamda**2*np.transpose(qef).dot(qef)) # !!!!!!!!!!!!!! A refaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DeltaL = DeltaL0
        lamda0 = lamda
        #Update the truss:
        self.truss.update()
        self.archive(lamda, current_iteration)
        self.display(step, lamda) 
        
        #####################
        ## PREDICTOR PHASE ##
        #####################
        while lamda < self.lamdaMax and not self.stopit:
            print "-- NormalPlaneArcLengthAlgorithm:: The predictor at step ",step,"     \r",
            step+=1
            #Apply loads to the truss:
            self.truss.applyNodalLoads(lamda)
            
            current_iteration = 0
            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            #Get the tangent stiffness matrix:
            Kt = self.truss.buildKt(self.computeTangentMethod)
            #Compute the out of balance forces:
            g = self.truss.getOOBF()
            #Compute the global loads vector:
            qef = self.truss.get_qef() 
            #Try to invert the tangent stiffness matrix:
            try:
                Kt_inv = np.linalg.inv(Kt)
            except numpy.linalg.LinAlgError:
                print "Kt is not invertible"
                print(Kt)
                sys.exit()
            else:
                deltaP  = - Kt_inv.dot(g)
                deltaPt =   Kt_inv.dot(qef)
                
                Ct1 = np.transpose(deltaP).dot(deltaP)*np.transpose(deltaPt).dot(deltaPt)+ (DeltaL**2-np.transpose(deltaP).dot(deltaP))*(np.transpose(deltaPt).dot(deltaPt)+self.psi**2*np.transpose(qef).dot(qef))
                if Ct1 < 0:
                    print "Ct1 is negative. Can't take its square root."
                    sys.exit()
                else:
                    Ct1 = math.sqrt(Ct1)
                #Verify that the matrix is positive definite:
                if np.all(np.linalg.eigvals(Kt) > 0): 
                    dlamda = (-np.transpose(deltaP).dot(deltaPt) + Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                else:
                    dlamda = (-np.transpose(deltaP).dot(deltaPt) - Ct1) / (np.transpose(deltaPt).dot(deltaPt) + self.psi**2 * np.transpose(qef).dot(qef))
                dLamdaPredict = dlamda
                du = deltaP + dlamda*deltaPt
                DuPredict = du
                Delta_u+=du
                lamda = lamda0 + dlamda
                #Update the positions:
                self.truss.incrementPositions(Delta_u)
                #Compute the error:
                error = self.computeError(g, Kt)
                
                #####################
                ## CORRECTOR PHASE ##
                #####################                
                while error > self.toll and current_iteration < self.nItMax:
                    #print "lamda = ", lamda
                    #Apply loads:
                    self.truss.applyNodalLoads(lamda)
                    #Build Kt:
                    Kt = self.truss.buildKt(self.computeTangentMethod)
                    #Get the out of balance forces:
                    g0 = self.truss.getOOBF()
                    #Get the global loads vector:
                    qef = self.truss.get_qef()
                    #Try to invert Kt:
                    try:
                        Kt_inv = np.linalg.inv(Kt)
                    except numpy.linalg.LinAlgError:
                        print "Kt is not invertible"
                        print(Kt)
                        sys.exit()
                    else:
                        deltaP = - Kt_inv.dot(g0) 
                        deltaPt =  Kt_inv.dot(qef)                    
                        #Solve equation (18) by first finding both coefficients:
                        a1 = dLamdaPredict*self.psi**2*np.transpose(qef).dot(qef)
                        a1 = a1 + np.transpose(deltaPt).dot(DuPredict)
                        a2 = np.transpose(DuPredict).dot(deltaP)

                        lamda1 = -a2/a1
                        dlamda += lamda1

                        du = deltaP + lamda1*deltaPt
                        #Update:
                        Delta_u += du
                        lamda = lamda0 + dlamda 
                        self.truss.incrementPositions(Delta_u)
                        error = self.computeError(g0, lamda)
                        current_iteration += 1
                        
                       
                if current_iteration != 0:     
                    DeltaL = DeltaL0*math.sqrt(float(self.Id)/current_iteration)

                else:
                    print "Error : should not end up here !"
                self.truss.update()
                self.archive(lamda[0,0], current_iteration)
                self.display(step, lamda)
                lamda0 = lamda
                
                
        print ""
        return 1