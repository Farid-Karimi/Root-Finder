# Rootfinding Application

Welcome to the Rootfinding Application! This bonus project, undertaken during my basic numerical analysis course, aims to offer a user-friendly interface for implementing diverse numerical analysis methods to find roots of mathematical functions. The application uses the Qt framework for the user interface, while the root finding methods are implemented using Python.

<details>
  <summary><h3>click here to see the images</h3></summary>

### Input page:
![ui](https://github.com/Farid-Karimi/Root-Finder/assets/118434072/4360d490-bfb3-4b21-8e54-a43e7e2a78c1)
### Result page:
![res](https://github.com/Farid-Karimi/Root-Finder/assets/118434072/3bef6e9b-f76b-443b-ae7e-198f52b17449)
### Plot interface:
![plot](https://github.com/Farid-Karimi/Root-Finder/assets/118434072/70beaf58-cf6e-4f6a-b7ab-5a06fd170d11)

</details>

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
3. Launch the application using your preferred Python IDE or by running the main script (pycharm is advised).
4. Start using the application by following the instructions provided above.

## Conclusion

The Rootfinding Application simplifies the process of finding roots of mathematical functions by providing a user-friendly interface and implementing various numerical analysis methods.
For future enhancements, the application could incorporate modified versions of the already implemented methods i.e. i struggled to implement the improved version of newton method for double roots so that could be a great starting point to improve the application.
