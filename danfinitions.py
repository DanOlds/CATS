# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 10:44:49 2016

@author: d5o
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def top_down_plot(data,**kwargs):    
    """plots 2D data from a pandas dataframe.  Requires that both index and column values of data frame be numeric!
    Optional arguments are colormap minimum and maximum (cmin, cmax), colormap to use (cmap), x-min and x-max 
    (xmin, xmax), y-value min and max (ymin, ymax),  if you want gouraud-style blurring (blur = True, default False),
    and if you want a colorbar plotted (colorbar = True, default False)"""
    import matplotlib.pyplot as plt

    X,Y = np.meshgrid(data.index,data.columns)
    X=X.T
    Y=Y.T

    cmin = min(data.min())
    cmax = max(data.max())
    use_cmap = 'viridis'
    xmin = min(data.index)
    xmax = max(data.index)
    ymin = min(data.columns)
    ymax = max(data.columns)
    shading_choice = 'None'
    use_colorbar = False
    
    if kwargs is not None:
        for key, value in kwargs.items():
            #print ("%s set to %s" %(key,value))
            if key == 'xmin':
                xmin = value
            if key == 'xmax':
                xmax = value
            if key == 'ymax':
                ymax = value
            if key == 'ymin':
                ymin = value
            if key == 'cmin':
                cmin = value
            if key == 'cmax':
                cmax = value
            if key == 'cmap':
                use_cmap = value
            if key == 'blur':
                if value == True:
                    shading_choice = 'gouraud'
                else:
                    shading_choice = 'None'
            if key == 'colorbar':
                if value == True:
                    use_colorbar = True
    plt.plot
    plt.pcolormesh(X,Y,data,cmap=use_cmap,vmin = cmin, vmax = cmax, shading = shading_choice)
    plt.axis([xmin,xmax,ymin,ymax])
    if use_colorbar:
        plt.colorbar()
    plt.show()

def is_number(arg1):
    if len(arg1.strip())>0:
        try:
            float(arg1)
            return True
        except ValueError:
            return False

def calc_rw(indata,infit,return_crw = False):
    fit = np.array(infit)
    data = np.array(indata)
    
    if return_crw == False:
        top = sum((data-fit)**2.0)
        bottom = sum(data**2.0)
        if bottom != 0.0:
            return ((top/bottom)**0.5)
        else:
            return 0.0
    else:
        vals = np.zeros(len(data))
        tdata_sum = 0.0
        tdiff_sum = 0.0
        for i in range(len(data)):
            tdata_sum += (data[i])**2
            tdiff_sum += (data[i] - fit[i])**2
            if tdata_sum != 0.0:
                vals[i] = tdiff_sum/tdata_sum
            else:
                vals[i] = 0.0
        vals = vals**(0.5)
        return vals
    
            
def calc_residual(y1,y2,return_sum = True,return_abs = True):
    res = np.zeros(len(y1))
    if  return_abs == False:
        res = (y1 - y2)
    if return_abs == True:
        res = abs(y1 - y2)
        
    if return_sum:
        return res.sum()
    else:
        return res
            
def read_in_all_runs_smart_list(arg1):
    xdata = []
    ydata = []    
    number_of_runs = len(arg1)

    print "reading in this many runs: "+str(number_of_runs)

    for j in range(number_of_runs):
        run_name = arg1[j]
        print "file reading in: "+str(run_name)
        with open(run_name, "r") as data:
            read_file = data.readlines()
            #setup arrays
            xvals = []
            yvals = []
            
            #determine own junk length
            for i in range(len(read_file)):
                if len(read_file[i].strip()) > 0:
                    if is_number(read_file[i].split()[0]):
                        #print "I want to start here "+str(i)
                        is_good = True
                        #now check consecutive ones
                        for k in range(i,i+10):
                            if(len(read_file[k].strip())) > 0:
                                if is_number(read_file[k].split()[0]):
                                    pass#still good
                                else:
                                    is_good = False
                            else:
                                is_good = False
                        if is_good:
                            #print "found a start "+str(i)
                            junk_length = i - 1
                            break
                            
            data_length = len(read_file)-junk_length

            
            #using found junk-length, assign data to xvals and yvals
            for i in range(junk_length+1,len(read_file)):
                xwouldbe = read_file[i].split()[0]
                if(is_number(xwouldbe)):                
                    xvals.append(float(read_file[i].split()[0]))
                    yvals.append(float(read_file[i].split()[1]))

            xdata.append(xvals)
            ydata.append(yvals)
            
    return xdata, ydata

def calculate_run_to_run_differences(arg1):
#calculate the differences between runs
    run_to_run_difference = []
    ydata = arg1
    number_of_runs = len(ydata)
    for j in range(number_of_runs-1):
        this_difference = np.subtract(ydata[j+1],ydata[j])
        run_to_run_difference.append(this_difference)
    return run_to_run_difference

def calculate_area_under_curves(arg1):    
#calculate the total area under the curve for each dataset
    ydata = arg1
    number_of_runs = len(ydata)
    sum_area_under_curve = []
    for j in range(number_of_runs):
        this_sum = []
        sum_is = 0.0
        for i in range(len(ydata[j])):
            sum_is += ydata[j][i]
            this_sum.append(sum_is)
        sum_area_under_curve.append(this_sum)
    return sum_area_under_curve

def save_all_differences(arg1,arg2,arg3):
    xdata = arg1
    ydata = arg2
    run_to_run_difference = arg3
    number_of_runs = len(ydata)
#save this file-to-file differences in a file
    outfile_name = "a_differences_run_to_run.dat"
    outfile2_name = "a_weighted_diff_run_to_run.dat"
    outfile = open(outfile_name,'w')
    outfile2 = open(outfile2_name,'w')
    offset = 0
    for j in range(number_of_runs-2):
        for i in range(len(xdata[j])):
            val = abs(run_to_run_difference[j][i])/(abs(ydata[j][i])+abs(ydata[j+1][i]))
            outfile.write(str(xdata[j][i])+" "+str(run_to_run_difference[j][i])+"\n")
            outfile2.write(str(xdata[j][i])+" "+str(val)+"\n")
        offset += 0.2
        outfile.write("\n")
        outfile2.write("\n")
    outfile.close()
    outfile2.close()
    
def save_all_data_in_one_file(arg1,arg2):
    xdata = arg1
    ydata = arg2
#save all data to a single file
    outfile_name = "a_all_data.dat"
    outfile2 = open("a_3dall_data.dat",'w')
    outfile = open(outfile_name,'w')
    for j in range(len(xdata)):
        for i in range(len(xdata[j])):
            outfile.write(str(xdata[j][i])+" "+str(ydata[j][i])+"\n")
            outfile2.write(str(j)+" "+str(xdata[j][i])+" "+str(ydata[j][i])+"\n")
        outfile2.write("\n")
        outfile.write("\n")
    outfile.close()
    outfile2.close()
    
def save_area_under_curves(arg1,arg2):
    xdata = arg1
    sum_area_under_curve = arg2
    outfile_name = "a_area_under_curve_sum.dat"
    outfile = open(outfile_name,'w')
    for j in range(len(xdata)):
        for i in range(len(xdata[j])):
            outfile.write(str(xdata[j][i])+" "+str(sum_area_under_curve[j][i])+"\n")
        outfile.write("\n")
    outfile.close()   


def calculate_diff_of_sum_vs_r(arg1):
#calculate the difference sum of targets as a function of r
    targetdata = arg1
    target_sum = np.copy(targetdata)
    for j in range(len(targetdata)):
        this_sum = 0.0
        for i in range(len(targetdata[0])):
            this_sum += targetdata[j][i]
            target_sum[j][i] = this_sum
    return target_sum

def calculate_scores_vs_r(arg1,arg2,arg3):
    number_of_runs = arg1
    targetdata = arg2
    ydata = arg3
#now calculate the scores of the datasets as a function of r (no scaling)
    fit_to_target = np.zeros([number_of_runs,len(targetdata)])
    fit_to_target_sum = np.zeros([number_of_runs,len(targetdata),len(ydata[0])])
    for k in range(len(targetdata)):
        for j in range(number_of_runs):
            this_goodness = 0.0
            for i in range(len(ydata[0])):
                this_goodness += abs((targetdata[k][i])-ydata[j][i])
                fit_to_target_sum[j][k][i] = this_goodness
            fit_to_target[j][k] = this_goodness
    return fit_to_target, fit_to_target_sum

def save_scores_vs_r(arg1):
    fit_to_target = arg1
    outfile_name = "a_scores_to_targets.dat"
    outfile = open(outfile_name,'w')
    for k in range(len(fit_to_target[0])): # 0 to 4
        for j in range(len(fit_to_target)): # 0 to 228
            outfile.write(str(j)+" "+str(fit_to_target[j][k])+"\n")
        outfile.write("\n")
    outfile.close()

def save_scores_to_target_sum(arg1,arg2):   
    xdata = arg1
    fit_to_target_sum = arg2
    number_of_runs = len(xdata)
    for k in range(len(fit_to_target_sum[0])):
        outfile_name = "a_scores_to_target_"+str(k)+"_sum_vs_r.dat"
        outfile = open(outfile_name,'w')
        for j in range(number_of_runs):
            for i in range(len(xdata[0])):
                outfile.write(str(xdata[0][i])+" "+str(fit_to_target_sum[j][k][i])+"\n")
            outfile.write("\n")
        outfile.close()
    
def make_average(arg1,arg2):
    ydata = arg1
    avg_list = arg2
    avg_data = np.zeros(len(ydata[0]))
    for j in avg_list:    
        for i in range(len(ydata[0])):
            avg_data[i] += ydata[j][i]
    avg_data[:] = avg_data[:] / (len(avg_list))
    return avg_data
    
    
def make_std_dev(arg1,arg2,arg3):
    ydata = arg1
    avg_list = arg2
    avg_data = arg3
    std_dev_list = np.zeros(len(ydata[0]))    
    for i in avg_list:
        for j in range(len(ydata[0])):
            std_dev_list[j] += (avg_data[j]-ydata[i][j])**2
    std_dev_list[:] = std_dev_list[:] / len(avg_list)
    std_dev_list[:] = std_dev_list[:]**(0.5)
    return std_dev_list
        

def try_all_the_things(mom,dad,df,num_pts=101,confidence = 0.05,return_error_bars = True,return_full_res_map = False):
        
    def two_parent_combo(mom,mom_frac,dad,dad_frac):
        return mom*mom_frac + dad*dad_frac
    df_phi_score = pd.DataFrame(index=np.linspace(0,1,num_pts),columns=df.columns,data=0)

    
    best_phi_list = []
    best_res_list = []

    if confidence > 1.0:
        confidence *= .01
    
    for cols in df_phi_score.columns:
        
        best_phi = 0
        best_res = 1e15
        
        test_data = df.loc[:,cols]
        for rows in df_phi_score.index:
            fit_frac = rows
            #this_res = calc_residual(two_parent_combo(mom,fit_frac,dad,1-fit_frac),test_data)
            this_res = calc_rw(two_parent_combo(mom,fit_frac,dad,1-fit_frac),test_data)

            if this_res < best_res:
                best_res = this_res
                best_phi = fit_frac

            df_phi_score.loc[rows,cols] = this_res
        
        best_phi_list.append(best_phi)
        best_res_list.append(best_res)
        
    best_phi_list = np.array(best_phi_list)
    best_res_list = np.array(best_res_list)
    
    wbp_low_list = []
    wbp_high_list = []
    
    if return_error_bars:
        for i in range(len(df_phi_score.columns)):
            cols = df.columns[i]
            worst_best_res = (1.0+confidence) * best_res_list[i]

            wbp_low = 0
            wbp_high = 1

            for j in range(0,len(df_phi_score.index),1):
                rows = df_phi_score.index[j]
                if df_phi_score.loc[rows,cols] <= worst_best_res:
                    wbp_low = rows
                    break

            for j in range(len(df_phi_score.index)-1,-1,-1):
                rows = df_phi_score.index[j]
                if df_phi_score.loc[rows,cols] <= worst_best_res:
                    wbp_high = rows
                    break

            wbp_high_list.append(wbp_high)
            wbp_low_list.append(wbp_low)



        wbp_high_list = np.array(wbp_high_list)
        wbp_low_list = np.array(wbp_low_list)

        if return_full_res_map == False:
            return best_phi_list, best_res_list, wbp_high_list, wbp_low_list
        else : #return the full res_map
            return best_phi_list, best_res_list, wbp_high_list, wbp_low_list,df_phi_score
    
    else : #don't bother with error bars
        if return_full_res_map == False:
            return best_phi_list, best_res_list
        else : #return the full res_map
            return best_phi_list,best_res,df_phi_score