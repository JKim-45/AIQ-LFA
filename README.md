## AIQ-LFA
R code for Early Quantitative Detection of Pancreatic Ductal Adenocarcinoma via an AI- and Quantum Dot Nanolens-powered Lateral Flow Assay: AIQ-LFA

## Paper Title
Early Quantitative Detection of Pancreatic Ductal Adenocarcinoma via an AI- and Quantum Dot Nanolens-powered Lateral Flow Assay: AIQ-LFA

## Authors
Han-Joo Bae, Minsup Shin, Hyunjoo Lee, Jun-Sik Chu, Kwanghee Yoo, Sohyeon Jang, Yuna Youn, Jaehyun An, Jin-Hyeok Hwang, Jaehi Kim,* Jong-chan Lee,* Luke P. Lee,* and Bong-Hyun Jun*


How to use these codes:  

Code for statistical analysis and machine learning (ML code_for submission.R)

## Raw data is in the 'raw datas.xlsx' file.

	stats 0 means sample from healthy control.
	stats 1 means sample from PDAC patient (row 52 to 91 : early-stage PDAC, row 92 to 151 : late-stage PDAC)
	EC means concentration of CA19-9 measured by using ELISA (AI-ELISA).
	QC means concentration of CA19-9 measured by using AIQ-LFA (AIQ-LFA).

## Before running, you can revise the following parameters by editing params.csv file; 

	ran_seed - random seed for classification of training set / validation set (in this manuscript, 42)
	cutoff - 0 : no dichotomization of CA19-9 value, 1 : dichotomization of CA19-9 value  
	kind_of_system - 0 : AI-ELISA, 1 : AIQ-LFA  
	n_of_boot - the number of bootstrapping (in this manuscript, 1000)
	train_ratio - the ratio of training set among whole sample (0 to 1 by 0.1 - in this manuscript, 0.7)
	FPR_rate - Using when finding the cutoff of probability (0 to 1, FPR < FPR_rate)
	TPR_rate - Using when finding the cutoff of probability (0 to 1, TPR > TPR_rate)

	Invalid value of parameters may cause error, please input valid parameters.

## Prior to run the code, you should install following packages; readxl, pROC, openxlsx, caret, and dplyr.

## Run the code after setting the parameters.

## When you run the code, then new folder will be created and data files will be generated in that folder;  

	HE = Healthy control versus Early PDAC, all = all samples  
	HL = Healthy control versus Early PDAC, but not used in this manuscript (just for supporting information)  

	AUC_test_~~.xlsx : AUC values with validation set  
	AUC_train_~~.xlsx : AUC values with training set  
	CI range.xlsx : 95% CI of AUC in each case  
	Cutoffs (FPR) at TPR more than TPR_rate.xlsx : Cutoffs at TPR > TPR_rate  
	Specificities at TPR more than TPR_rate.xlsx : Specificities at TPR > TPR_rate  
	Cutoffs (TPR) at FPR less than FPR_rate.xlsx : Cutoffs at FPR < FPR_rate  
	Sensitivities at FPR less than FPR_rate.xlsx : Sensitivities at FPR < FPR_rate  
	confusion_matrix_fpr_~~.xlsx : confusion matrix of validation set with cutoff from FPR < FPR_rate  
	confusion_matrix_tpr_~~.xlsx : confusion matrix of validation set with cutoff from TPR > TPR_rate  
	P values (Status from P).xlsx : P value from DeLong test of ROC curves from training set and validation set, in each repetition  
	Delong.xlsx : DeLong test of ROC curves from training set and validation set, and average P value  
	ROC_test_~~.xlsx : representative ROC curve from validation set  
	ROC_train_~~.xlsx : representative ROC curve from training set  
	test_for_~~.xlsx : validation set  
	train_for_~~.xlsx : training set  
	Parameters.xlsx : the parameter settings 


