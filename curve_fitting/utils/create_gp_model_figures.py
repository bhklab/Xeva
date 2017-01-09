"""
Quick script for creating GP figures for PDX data from R to
Pandas dataframes, using GPy and matplotlib.

Want a figure of the fit for each replica within each model.

So for each model and for each model's replicate, will want
three charts:

1. Raw data
2. Gaussian process fit onto control data with GP parameters
3. Gaussian process fit onto case data with GP parameters

"""

# R and Dataframe deps
import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas as pd

# Plotting deps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# GPy deps
import GPy
GPy.plotting.change_plotting_library('matplotlib')

# Logging
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def plot_original_data(model, growth_data, collection_days, replicate,
                       control_or_case):
    """
    Plots and saves the original timeseries data for the
    replicate.

    :param model:
    :param growth_data:
    :param collection_days:
    :param replicate:
    :param control_or_case:
    :return:
    """

    df_PDXNNN_GrowthData = pd.DataFrame(growth_data).dropna()
    df_PDXNNN_GrowthData["Days"] = collection_days[:len(df_PDXNNN_GrowthData)]

    plot = df_PDXNNN_GrowthData.plot(y=replicate, x='Days', kind='scatter')
    plot.set_title("Original Data of " + model + " " + control_or_case + ", " +
                     "Replicate " + replicate)
    plot.set_xlim(-10, 150)
    plot.set_ylim(0, 2200)
    plot.set_xlabel('Time (days)')
    plot.set_ylabel('Tumor size (mm^3)')

    figure = plot.get_figure()
    return figure


def plot_fitted_data(model, growth_data, collection_days, replicate,
                     control_or_case):
    """
    Plots and saves the Gaussian process-fitted timeseries
    data for replicate.

    :param model:
    :param growth_data:
    :param collection_days:
    :param replicate:
    :param control_or_case:
    :return:
    """

    df_PDXNNN_GrowthData = pd.DataFrame(growth_data).dropna()
    df_PDXNNN_GrowthData["Days"] = collection_days[:len(df_PDXNNN_GrowthData)]

    # kernel setting and model optimization
    kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=5.)
    X = np.asarray([[day] for day in df_PDXNNN_GrowthData['Days']])
    Y = np.asarray([[size] for size in df_PDXNNN_GrowthData[replicate]])
    gp_model = GPy.models.GPRegression(X, Y, kernel)
    gp_model.optimize_restarts(num_restarts=10, verbose=False)

    # plotting
    tau = 0.05
    grid_params = [-25, max(X) + 25, tau]
    extent = np.array(grid_params[:2])

    plot = gp_model.plot(plot_limits=extent)
    plot.set_title('Gaussian Process fit of ' + model + ' ' + control_or_case +
                     ',' + ' Replicate ' + replicate)
    plot.set_xlim(extent)
    plot.set_ylim(0, 2200)
    plot.set_xlabel('Time (days)')
    plot.set_ylabel('Tumor size (mm^3)')
    plot.text(-10, 2100, 'RBF Variance: ' + str(gp_model.param_array[0]))
    plot.text(-10, 2000, 'RBF Lengthscale: ' + str(gp_model.param_array[1]))
    plot.text(-10, 1900, 'Gaussian noise variance: ' + str(gp_model.param_array[2]))

    figure = plot.get_figure()
    return figure

def process_all_control_replicates(model, original_control_save_to_path,
                                   fitted_control_save_to_path):
    """
    Main processing function for one model. Gets all control replicates
    for model, plots original and fits a GP.

    :param model:
    :param original_control_save_to_path:
    :param fitted_control_save_to_path:
    :return:
    """

    # load collection days and drug start day
    PDXNNN_CollectionDays = com.load_data(model + '_CollectionDays')
    PDXNNN_DrugStartDay = com.load_data(model + '_DrugStartDay')

    # control first
    PDXNNN_ControlGrowth = com.load_data(model + '_ControlGrowth')
    control_replicates = PDXNNN_ControlGrowth.columns

    failed_replicates = []

    logger.info("Processing control replicates.")
    for replicate in control_replicates:

        logger.info("Plotting original control data...")
        original_data_control_fig = plot_original_data(model, PDXNNN_ControlGrowth,
                                                       PDXNNN_CollectionDays,
                                                       replicate, 'Control')
        original_control_save_to_path.savefig(original_data_control_fig)
        plt.close(original_data_control_fig)

        logger.info("Plotting fitted control data...")
        try:
            fitted_data_control_fig = plot_fitted_data(model, PDXNNN_ControlGrowth,
                                                       PDXNNN_CollectionDays,
                                                       replicate, 'Control')
            fitted_control_save_to_path.savefig(fitted_data_control_fig)
            plt.close(fitted_data_control_fig)
        except Exception as e:
            logger.info(model + " " + replicate + " failed fit process.")
            logger.info(e)
            failed_replicates.append((model, replicate))

    return failed_replicates

def process_all_case_replicates(model, original_case_save_to_path, fitted_case_save_to_path):
    """
    Main processing function for one model. Gets all case replicates
    for model, plots original and fits a GP.

    :param model:
    :param original_case_save_to_path:
    :param fitted_case_save_to_path:
    :return:
    """

    # load collection days and drug start day
    PDXNNN_CollectionDays = com.load_data(model + '_CollectionDays')
    PDXNNN_DrugStartDay = com.load_data(model + '_DrugStartDay')

    # cases
    PDXNNN_CaseGrowth = com.load_data(model + '_CaseGrowth')
    case_replicates = PDXNNN_CaseGrowth.columns

    failed_replicates = []

    logger.info("Processing case replicates.")
    for replicate in case_replicates:
        logger.info("Plotting original case data...")
        original_data_case_fig = plot_original_data(model, PDXNNN_CaseGrowth,
                                                    PDXNNN_CollectionDays,
                                                    replicate, 'Case')
        original_case_save_to_path.savefig(original_data_case_fig)
        plt.close(original_data_case_fig)

        logger.info("Plotting fitted case data...")
        try:
            fitted_data_case_fig = plot_fitted_data(model, PDXNNN_CaseGrowth,
                                                    PDXNNN_CollectionDays,
                                                    replicate, 'Case')
            fitted_case_save_to_path.savefig()
            plt.close(fitted_data_case_fig)
        except:
            logger.info(model + " " + replicate + " failed fit process.")
            failed_case_replicates.append((model, replicate))

    return failed_replicates


if __name__ == '__main__':

    # load data
    robjects.r['load']("data/readyForPandasKrasp53.RData")

    # declare all plots to be made
    to_process_file = open('utils/KRASp53ToProcess.txt', 'r')
    original_control_save_to = PdfPages('figures/KRASp53PlottedOriginalControls.pdf')
    fitted_control_save_to = PdfPages('figures/KRASp53PlottedGPControls.pdf')
    original_case_save_to = PdfPages('figures/KRASp53PlottedOriginalCase.pdf')
    fitted_case_save_to = PdfPages('figures/KRASp53PlottedGPCase.pdf')

    # GPy may sometimes raise a GPy error - this
    # holds all of the model-replicate names
    failed_control_replicates = []
    failed_case_replicates = []

    for model in to_process_file:
        failed_control_replicates += \
            process_all_control_replicates(model.strip(),
                                           original_control_save_to,
                                           fitted_control_save_to)
        failed_case_replicates += \
            process_all_case_replicates(model.strip(),
                                        original_case_save_to,
                                        fitted_case_save_to)
    logger.info("Failed control replicates:")
    print(failed_control_replicates)
    logger.info("Failed case replicates:")
    print(failed_case_replicates)

    original_control_save_to.close()
    fitted_control_save_to.close()
    original_case_save_to.close()
    fitted_case_save_to.close()