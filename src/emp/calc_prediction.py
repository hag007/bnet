import sys
sys.path.insert(0, '../..')
import json
import argparse
import os
import src.constants as constants
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn import svm

def svm_rbf_default(tuned_parameters):
    return GridSearchCV(svm.SVC(probability=True), param_grid=tuned_parameters, return_train_score=True)

def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', help='/path/to/network_file', default=constants.config_json["network_file"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default=constants.config_json["additional_args"])
    parser.add_argument('--phenotype_args', help="phenotype_args", dest='phenotype_args',
                        default=constants.config_json["phenotype_args"])
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    network_file = args.network_file
    true_solutions_folder = args.true_solutions_folder
    phenotype_args = args.phenotype_args # json.loads(args.phenotype_args)
    additional_args = args.additional_args # json.loads(args.additional_args)

    calc_prediction(dataset_file, algo, network_file, true_solutions_folder, additional_args, phenotype_args)


def calc_prediction(dataset_file, algo, network_file, true_solutions_folder, additional_args, phenotype_args):

    phenotype_args = json.loads(phenotype_args)
    additional_args = json.loads(additional_args)

    params_name = "_".join([str(additional_args[a]) for a in \
                                ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor", "activity_baseline"]])

    # read files
    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    network_name=os.path.splitext(os.path.split(network_file)[1])[0]
    unique_folder_name="{}_{}_{}_{}".format(dataset_name,network_name,algo,params_name)
    output_folder=os.path.join(true_solutions_folder, unique_folder_name)
    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    output_file = os.path.join(output_folder,"report.tsv")

    pcs_data=pd.read_csv(os.path.join(output_folder, "pcs.tsv"), sep='\t', index_col=0)
    if pcs_data.shape[0]==0:
        with open(output_file, 'w') as o:
            print(f'conf: {additional_args}\n{output_file}:\n no PCS here...', 'r')
            o.write("\t".join(["name", "N features", "RF accuracy","SVM accuracy","SVM AUPR","SVM AUROC","Null AUPR","Null AUROC"])+"\n")
            o.write("\t".join([unique_folder_name, "0", "0", "0", "0", "0", "0", "0"]))
            return

    phenotype = pd.read_csv(phenotype_args["path to TCGA file"], sep='\t')

    # parse data
    phenotype.index = phenotype.loc[:,constants.config_json["phenotype_index"]]
    phenotype_field = phenotype_args["name of variable"]
    phenotype=phenotype.loc[:,phenotype_field]
    phenotype_0=phenotype[phenotype.isin(phenotype_args["value1"])]
    phenotype_1 = phenotype[phenotype.isin(phenotype_args["value2"])]
    phenotype=pd.concat([phenotype_0,phenotype_1], axis=0)

    phenotype=phenotype.reindex(pcs_data.index).dropna(axis=0)
    pcs_data = pcs_data.reindex(phenotype.index).dropna(axis=0)
    # for a in [phenotype_args["value1"],phenotype_args["value2"]]:
    #     print(f'# of phenotype {a}: {(phenotype.isin(a)).sum()}')

    # prepare structures
    prediction_accuracies_rf=[]
    prediction_accuracies_svm = []
    svm_pr_aucs=[]
    svm_roc_aucs=[]
    svm_pr_aucs_nulls=[]
    svm_roc_aucs_nulls=[]

    rf_pr_aucs = []
    rf_roc_aucs = []
    rf_pr_aucs_nulls = []
    rf_roc_aucs_nulls = []

    # calc prediction n_iterations times

    labels = [([a in phenotype_args["value1"], a in phenotype_args["value2"]].index(True)) for a in phenotype]

    n_iterations=100
    for a in np.arange(n_iterations):
        # run RF
        X_train, X_test, y_train, y_test = train_test_split(pcs_data, labels, test_size = 0.33)
        clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf.fit(X_train, y_train)
        y_pred=clf.predict(X_test)
        probabilities = clf.predict_proba(X_test)
        probabilities=[a[1] for a in  probabilities]
        rf_pr_aucs.append(average_precision_score(y_test, probabilities))
        rf_roc_aucs.append(roc_auc_score(y_test, probabilities))

        # calc "null" metrics
        rf_pr_aucs_nulls.append(average_precision_score(y_test, [1 for a in y_test]))
        rf_roc_aucs_nulls.append(roc_auc_score(y_test, [1 for a in y_test]))

        prediction_accuracies_rf.append(float(np.sum([a==b for a,b in zip(y_test, y_pred)]))/len(y_test))  # roc_auc_score(labels_test, probabilities))

        # run SVM
        classification_method = "svm_rbf_default"
        tuning_parameters = {'C': [10], 'kernel': ['rbf']}
        thismodule = sys.modules[__name__]
        clf_method = getattr(thismodule, classification_method)(tuning_parameters)
        clf_method.fit(X_train, y_train)
        predicted_results = clf_method.predict(X_test)

        # calc metrics (AUROC, AUPR)
        probabilities = clf_method.decision_function(X_test)
        # probabilities = [cur[predicted_results[i]] for i, cur in enumerate(clf_method.predict_proba(X_test))] // RAnking alternative. Has some issues...
        precision, recall, _ = precision_recall_curve(y_test, probabilities)
        svm_pr_aucs.append(average_precision_score(y_test, probabilities))
        svm_roc_aucs.append(roc_auc_score(y_test, probabilities))

        # calc "null" metrics
        svm_pr_aucs_nulls.append(average_precision_score(y_test, [1 for a in y_test]))
        svm_roc_aucs_nulls.append(roc_auc_score(y_test, [1 for a in y_test]))

        prediction_accuracies_svm.append(float(np.sum([a == b for a, b in zip(y_test, predicted_results)])) / len(
            y_test))  # roc_auc_score(labels_test, probabilities))

    # print averaged metrics
    # print(f'RF accuracy: {np.mean(prediction_accuracies_rf)}')
    # print(f'SVM accuracy: {np.mean(prediction_accuracies_svm)}')
    # print(f'SVM average AUPR: {np.mean(svm_pr_aucs)} (AUPR null (all labels are 1): {np.mean(svm_pr_aucs_nulls)})')
    # print(f'SVM average AUROC: {np.mean(svm_roc_aucs)} (null score (all labels are 1) is : {np.mean(svm_roc_aucs_nulls)})')


    with open(output_file, 'w') as o:
        metrics = [prediction_accuracies_rf, rf_pr_aucs, rf_roc_aucs, rf_pr_aucs_nulls,
                   rf_roc_aucs_nulls, prediction_accuracies_svm, svm_pr_aucs, svm_roc_aucs, svm_pr_aucs_nulls,
                   svm_roc_aucs_nulls]
        o.write("\t".join(["name", "N features", "RF accuracy", "RF AUPR", "RF AUROC", "Null RF AUPR", "Null RF AUROC", "SVM accuracy", "SVM AUPR", "SVM AUROC", "Null SVM AUPR", "Null SVM AUROC"]) + "\n")
        o.write("\t".join([unique_folder_name, str(pcs_data.shape[1])] + [str(round(np.mean(m), 4)) for m in metrics]))

    with open(output_file, 'r') as r:
        print(f'conf: {additional_args}\npheno count:{np.sum(labels)},{len(labels)-np.sum(labels)}\n{output_file}:\n{r.read()}', 'r')




if __name__ == "__main__":
    main()
