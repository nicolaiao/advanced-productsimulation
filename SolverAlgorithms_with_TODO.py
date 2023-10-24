import math
import numpy as np
from CorotBeam_with_TODO import rot_matrix

def update_global_disp(v_hat, delta_v):
    for (dv_a, dw_a) in delta_v:
        rot = rot_matrix(dw_a)
        v_hat += dv_a #Will not work
        rot_a = rot + v_hat(rot_a) # will not work

    return v_hat # TODO: FIX THIS

def solveArchLength(model, archLength=0.02, max_steps=50, max_iter=30):
    """_summary_

    Args:
        model (_type_): _description_
        archLength (float, optional): _description_. Defaults to 0.02.
        max_steps (int, optional): _description_. Defaults to 50.
        max_iter (int, optional): _description_. Defaults to 30.
    
    Internals:
        uVec (np.array(num_dofs))

    """
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs) 
    res_Vec = np.zeros(num_dofs)
    Lambda = 0.0

    d_q_prev = np.zeros(num_dofs)

    for iStep in range(max_steps):

        #TODO: Implement this
        w_q0 = model.get_incremental_load(Lambda) # double check ??
        q = np.dot(model.get_K_sys(uVec), w_q0) # not correct uvec
        f = np.sqrt(1 + np.dot(np.transpose(w_q0), w_q0))
        if np.dot(np.transpose(w_q0), uVec) > 1:
            d_lambda = archLength/f
        else:
            d_lambda = -archLength/f

        uVec = d_lambda * w_q0 #v0 ??
        Lambda += d_lambda
        v_hat = update_global_disp(v_hat, delta_v) #What is v_hat and delta_v ????

        for iIter in range(max_iter):

            # TODO: Implement this
            #K(vˆ)wq = q
            #K(vˆ)wr = −r(vˆ, λ)
            #del_λ =  -(w_q0.T * w_r)/(1 + w_q0.T*w_q)
            res_Vec = model.get_residual(Lambda, uVec) # <-- Has to be this
            w_r = model.get_external_load(Lambda)
            tmp = np.dot(np.transpose(w_q0), w_r)
            d_lambda = - tmp / (1 + tmp)
            Lambda += d_lambda
            #node_disp #v_hat
            #res_Vec = model.get_residual(uVec, Lambda) #get_residual(load_factor, disp_sys)
            if (res_Vec.dot(res_Vec) < 1.0e-15): # until ||r(v, lambda)|| < epsilon
                break
            res_Vec = model.get_residual(uVec, Lambda) #Lambda has to be iterable
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                break

        model.append_solution(Lambda, uVec)
        print("Linear arc step {:}  load_factor= {:12.3e}".format(iStep, Lambda))

def solveNonlinLoadControl(model, load_steps=0.01, max_steps=100, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec   = np.zeros(num_dofs)
    d_uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep # Lambda: float, disp_sys: iterable, load_factor: float

        #TODO: Implement this #(Predictor step with arc length deltaS) ?
        for iIter in range(max_iter):
            # TODO: Implement this #(Corrector iterations) ?
            #K(vˆ)wq = q
            #K(vˆ)wr = −r(vˆ, λ)
            #del_λ =  -(w_q0.T * w_r)/(1 + w_q0.T*w_q)
            res_Vec = model.get_residual(Lambda, uVec) # <-- Has to be this
            w_q0 = model.get_incremental_load(Lambda) #This is w_q FIX
            w_r = model.get_external_load(Lambda)
            tmp = np.dot(np.transpose(w_q0), w_r)
            d_lambda = - tmp / (1 + tmp)
            Lambda += d_lambda
            #node_disp #v_hat
            #res_Vec = model.get_residual(uVec, Lambda) #get_residual(load_factor, disp_sys)
            if (res_Vec.dot(res_Vec) < 1.0e-15): # until ||r(v, lambda)|| < epsilon
                break




        model.append_solution(Lambda, uVec)
        print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))



def solveLinearSteps(model, load_steps=0.01, max_steps=100):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        q_Vec   = model.get_incremental_load(Lambda)

        K_mat = model.get_K_sys(uVec)

        d_q = np.linalg.solve(K_mat, q_Vec)

        uVec = d_q * Lambda

        model.append_solution(Lambda, uVec)
        print("Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))

