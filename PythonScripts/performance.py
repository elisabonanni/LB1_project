'''Calculate CONFUSION MATRIX and the PERFORMANCE providing as input a threshold'''

#!/usr/bin/env python3

import sys
import numpy as np

def get_data(predfile):
    '''Function that reads the file'''
    preds=[] #initialise a list
    f = open(predfile)
    for line in f:
        v=line.rstrip().split() #return identifiers, e.value and label
        v[1]=float(v[1]) #E-Value
        v[2]=int(v[2]) #Label
        preds.append(v)
    return preds

#CONFUSION MATRIX
def get_cm(preds,th=0.5): #take as input the preds list generated and threshold
    '''Function that out of the file calculate the Confusion Matrix'''
    cm=np.zeros((2,2)) #initialise a 2x2 matrix
    for pred in preds: #go line by line and assign the different cases (0,1) in the the cm
        p=0 #label 0, assigned by default
        r=pred[2]
        if pred[1]<=th:
            p=1 #if the value is lower than the threshodl assign 1
        cm[p][r]+=1 #predictions on rows and real values on columns
    return cm #return the confusion matrix TN, FN, FP, TP

#OVERALL ACCURACY
def get_accuracy (cm):
    '''Calculate the overall accuracy from the confusion matrix'''
    q2=float(cm[0][0]+cm[1][1])/np.sum(cm) #(TN+TP)/total
    return q2

#MATTHEW CORRELATION COEFFICIENT
def get_mcc(cm): #(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    '''Calculate Matthew Correlation Coefficient'''
    tp = cm[1][1]
    tn = cm[0][0]
    fp = cm[1][0]
    fn = cm[0][1]
    mcc = (tp*tn-fp*fn)/np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) #at the denominator all the possible combination between true and false
    return mcc

def get_f1(cm):
    tp = cm[1][1]
    tn = cm[0][0]
    fp = cm[1][0]
    fn = cm[0][1]
    prec= (tp/(tp+fp))
    rec= (tp/(tp+fn))
    f1 = (2*prec*rec)/(prec+rec)
    return f1

if __name__=='__main__':
  predfile=sys.argv[1]
  th=0.001
  if len(sys.argv)>2:
      th=float(sys.argv[2])
  preds=get_data(predfile)
  cm=get_cm(preds,th)
  q2=get_accuracy(cm)
  mcc=get_mcc(cm)
  f1=get_f1(cm)
  print('TH:',th,'Q2:',q2,'MCC:',mcc,'F1:',f1,'TP=',cm[1][1],'TN=',cm[0][0],'FP=',cm[1][0],'FN=',cm[0][1])
