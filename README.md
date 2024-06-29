# Rootfinding Application

Welcome to the Rootfinding Application! This project is designed to provide a user-friendly interface for implementing various numerical analysis methods to find roots of mathematical functions. The application utilizes the Qt framework for the user interface, while the root finding methods are implemented using Python.

## Methods Implemented

The Rootfinding Application incorporates the following numerical analysis methods:

1. Newton's Method
2. Bisection Method
3. Secant Method
4. Fixed Point Iteration Method
5. False Position Method

Each method is carefully implemented to accurately find the roots of the given function.

## How to Use

Follow the steps below to effectively utilize the application:

### Input Configuration

To start using the application, you need to provide the necessary input configuration. This includes the following:

1. **Function**: Enter the mathematical function for which you want to find the root. Make sure to input the function in a valid format.
2. **Error Tolerance**: Specify the desired level of accuracy for the root. This is optional, and if not provided, the default value of 0.0001 will be used.
3. **Number of Iterations**: Define the maximum number of iterations to perform. If not specified, the default value of 100 will be used.
4. **Plot function**: Users can input the interval for the function to be drawn in the `x_0` and `x_1` fields. Once the function and interval are provided, they can click the "Plot function" button to visualize an interactive plot of the function.

### Running the Application

Once you have entered the input configuration, simply click the "Calculate" button to initiate the root finding process. The application will execute the selected method and display the result on the user interface.

### Viewing the Result

After the root finding process is completed, the application will present the result on the interface. This includes the root value and any additional information, such as the number of iterations performed, how much time it took to find the root and the error in each iteration.

## Installation

To install and run the Rootfinding Application, please follow these steps:

1. Clone the repository to your local machine.
2. Install the required dependencies listed in the `requirements.txt` file.
3. Launch the application using your preferred Python IDE or by running the main script.
4. Start using the application by following the instructions provided above.

## Installation

To use the Rootfinding Application, follow these simple steps:

1. Clone the repository to your local machine.
2. Install the required dependencies listed in the `requirements.txt` file.
3. Launch the application using your preferred Python IDE or by running the main script.
4. Enjoy exploring and utilizing the various root finding methods provided!

## Conclusion

The Rootfinding Application simplifies the process of finding roots of mathematical functions by providing a user-friendly interface and implementing various numerical analysis methods.
For future enhancements, the application could incorporate modified versions of the already implemented methods i.e. i struggled to implement the improved version of newton method for double roots so that could be a great start to improve the application.
