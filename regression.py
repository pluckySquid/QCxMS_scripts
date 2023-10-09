import numpy as np
# from sklearn.linear_model import LinearRegression
# from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

# # Define the input data
# xData = np.array([13, 16, 20, 24, 28, 37])
# yData = np.array([13.4, 200, 121.8, 452.6, 1441.4, 3800])
# # xData = np.array([19.1647, 18.0189, 16.9550, 15.7683, 14.7044, 13.6269, 12.6040, 11.4309, 10.2987, 9.23465, 8.18440, 7.89789, 7.62498, 7.36571, 7.01106, 6.71094, 6.46548, 6.27436, 6.16543, 6.05569, 5.91904, 5.78247, 5.53661, 4.85425, 4.29468, 3.74888, 3.16206, 2.58882, 1.93371, 1.52426, 1.14211, 0.719035, 0.377708, 0.0226971, -0.223181, -0.537231, -0.878491, -1.27484, -1.45266, -1.57583, -1.61717])
# # yData = np.array([0.644557, 0.641059, 0.637555, 0.634059, 0.634135, 0.631825, 0.631899, 0.627209, 0.622516, 0.617818, 0.616103, 0.613736, 0.610175, 0.606613, 0.605445, 0.603676, 0.604887, 0.600127, 0.604909, 0.588207, 0.581056, 0.576292, 0.566761, 0.555472, 0.545367, 0.538842, 0.529336, 0.518635, 0.506747, 0.499018, 0.491885, 0.484754, 0.475230, 0.464514, 0.454387, 0.444861, 0.437128, 0.415076, 0.401363, 0.390034, 0.378698])

# import numpy, scipy, matplotlib
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
# from scipy.optimize import differential_evolution
# import warnings

# xData = np.array([19.1647, 18.0189, 16.9550, 15.7683, 14.7044, 13.6269, 12.6040, 11.4309, 10.2987, 9.23465, 8.18440, 7.89789, 7.62498, 7.36571, 7.01106, 6.71094, 6.46548, 6.27436, 6.16543, 6.05569, 5.91904, 5.78247, 5.53661, 4.85425, 4.29468, 3.74888, 3.16206, 2.58882, 1.93371, 1.52426, 1.14211, 0.719035, 0.377708, 0.0226971, -0.223181, -0.537231, -0.878491, -1.27484, -1.45266, -1.57583, -1.61717])
# yData = np.array([0.644557, 0.641059, 0.637555, 0.634059, 0.634135, 0.631825, 0.631899, 0.627209, 0.622516, 0.617818, 0.616103, 0.613736, 0.610175, 0.606613, 0.605445, 0.603676, 0.604887, 0.600127, 0.604909, 0.588207, 0.581056, 0.576292, 0.566761, 0.555472, 0.545367, 0.538842, 0.529336, 0.518635, 0.506747, 0.499018, 0.491885, 0.484754, 0.475230, 0.464514, 0.454387, 0.444861, 0.437128, 0.415076, 0.401363, 0.390034, 0.378698])

x_protonation_time = [0, 13, 16, 20, 24, 28, 37, 44, 54]
y_protonation_time = [0, 0.0076, 0.036, 0.0536, 0.0447, 0.0658, 0.137, 0.449, 0.39]
y_protonation_ram = [0, 0.0035, 0.062, 0.114, 0.1657, 0.333, 0.273, 0.647, 0.487]

x_qcxms_time = [0, 13, 16, 20, 24, 28, 37, 44, 54]
y_qcxms_time = [0, 0.075, 0.272, 0.282, 0.361, 0.5, 1.078, 1.278, 0.761]
y_qcxms_ram = [0, 0.4, 0.8, 1.02, 1.233, 1.843, 2.967, 3.9, 4.78]

x_qcxms_time = [0, 13, 16, 20, 24, 28, 37, 44, 54]
y_pqcxms_time = [0, 0.009, 0.18, 0.08, 0.247, 0.863, 1.985, 4.02, 3.816]
y_pqcxms_ram = [0, 0.05, 0.313, 0.246, 0.634, 1.72, 3.593, 5.621, 9.441]

x_pqcxms_time = [0, 13, 16, 20, 24, 28, 37, 44, 54]
y_total_time = [0, 0.3, 4.8, 2.4, 6.4, 24.8, 50.2, 104.3, 96.5]



xData = np.array([0, 13, 16, 20, 24, 28, 37,54,75,77])
yData = np.array([0, 13.4, 200, 121.8, 452.6, 1441.4, 3800, 20995.2, 64350, 72570])

xData = np.array(x_protonation_time)
yData = np.array(y_total_time)
#yData = np.array([6.152, 12.75, 16.63, 29.67, 44.83, 100])

z = np.poly1d(np.polyfit(xData,yData,2))
plt.plot(xData, yData)
x = np.arange(0,100)
y = z(x)
plt.plot(x, y)
plt.title('CPU-hours')
plt.xlabel('Atom (n)')
plt.ylabel('Time (h)')
plt.savefig('cpu.png')
plt.show()
# plt.plot(z)
# plt.show()
print(z)
# def func(x, a, b, Offset): # Sigmoid A With Offset from zunzun.com
#     return  1.0 / (1.0 + numpy.exp(-a * (x-b))) + Offset

# results: 
# protonation time: 0.0001765 x ^ 2 - 0.0009905 x - 0.007388
# qcxms time: 
# y_pqcxms_time: 0.001704 x ^ 2 - 0.005459 x - 0.2039
# y_total_time: 0.04182 x ^ 2 - 0.047 x - 5.641


















# # function for genetic algorithm to minimize (sum of squared error)
# def sumOfSquaredError(parameterTuple):
#     warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
#     val = func(xData, *parameterTuple)
#     return numpy.sum((yData - val) ** 2.0)


# def generate_Initial_Parameters():
#     # min and max used for bounds
#     maxX = max(xData)
#     minX = min(xData)
#     maxY = max(yData)
#     minY = min(yData)

#     parameterBounds = []
#     parameterBounds.append([minX, maxX]) # search bounds for a
#     parameterBounds.append([minX, maxX]) # search bounds for b
#     parameterBounds.append([0.0, maxY]) # search bounds for Offset

#     # "seed" the numpy random number generator for repeatable results
#     result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
#     return result.x

# # generate initial parameter values
# geneticParameters = generate_Initial_Parameters()

# # curve fit the test data
# fittedParameters, pcov = curve_fit(func, xData, yData, geneticParameters)

# print('Parameters', fittedParameters)

# modelPredictions = func(xData, *fittedParameters) 

# absError = modelPredictions - yData

# SE = numpy.square(absError) # squared errors
# MSE = numpy.mean(SE) # mean squared errors
# RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
# Rsquared = 1.0 - (numpy.var(absError) / numpy.var(yData))
# print('RMSE:', RMSE)
# print('R-squared:', Rsquared)



# ##########################################################
# # graphics output section
# def ModelAndScatterPlot(graphWidth, graphHeight):
#     f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
#     axes = f.add_subplot(111)

#     # first the raw data as a scatter plot
#     axes.plot(xData, yData,  'D')

#     # create data for the fitted equation plot
#     xModel = numpy.linspace(min(xData), max(xData))
#     yModel = func(xModel, *fittedParameters)

#     # now the model as a line plot 
#     axes.plot(xModel, yModel)

#     axes.set_xlabel('X Data') # X axis data label
#     axes.set_ylabel('Y Data') # Y axis data label

#     plt.show()
#     plt.close('all') # clean up after using pyplot

# graphWidth = 800
# graphHeight = 600
# ModelAndScatterPlot(graphWidth, graphHeight)

