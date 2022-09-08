# -*- coding: utf-8 -*-
"""
PHYS20161 Introduction To Programming For Physicists
Final Assignment: 'Doppler Spectroscopy'

This code takes in the data of the doppler shifts obtained from the star
during a period of several years, converts it to the velocity of the star
using a formula for Doppler's shift and fits a predicted curve to it by
minimising the chi squared.

The relative velocity of the Earth and the star is already subtracted from
the data.

A calculation is also performed to find the orbital distance and the mass of
the exoplanet, causing the observed variation in the velocity of the star, and
it's uncertainty.

Plots representing the data and the fit are also created to aid the analysis.

Ksenija Kovalenka ID:10485506
Date: 16/12/2020
"""

# libraries and predefined functions
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc
from scipy.optimize import fmin

# constants
FILE_NAMES = ['doppler_data_1.csv',
              'doppler_data_2.csv']
REST_WAVELENGHT = 656.281 # nm
STAR_MASS = 2.78 * 1.98847e+30 # kg
JOVIAN_MASS = 1.89813e+27 # kg


# functions
def read_and_combine(file_names, sort=True):
    """
    Combines two data sets. Outputs combined 2D arrays of floats.
    Optionally, sorts the data in ascending order by the first column.
    Parameters
    ----------
    file_names : array of strings
    sort : bool
    Returns
    -------
    combined_data : 2D array of floats
    file_found : bool

    """
    combined_data = np.empty((0, 3))

    for file in file_names:
        try:
            data_set = np.genfromtxt(file, delimiter=',', comments='%')
            file_found = True

        except OSError:
            print('File {} not found. PLease check the file names and put them'
                  ' in the same directory as the code.'.format(file))
            file_found = False
            break

        except ValueError:
            print('All lines must have the same number of elements (3 expected).'
                  ' Please check your data file.')
            file_found = False
            break

        try:
            combined_data = np.vstack((combined_data, data_set))

        except ValueError:
            print('Three columns of data expected. Please check your data file')
            file_found = False
            break

    if sort:
        combined_data = combined_data[combined_data[:, 0].argsort()]

    return combined_data, file_found


def validate(data_input):
    """
    1. Removes non-usable data, sush as Nans or infinities.
    2. Removes data with zero uncertainties, as such a measurement is not
       physically possible.
    3. Removes data too far from its neighbours, because values are expected
       to vary smoothly with time.
    Parameters
    ----------
    data_input : ND array of floats

    Returns
    -------
    good_data : ND array of floats

    """
    # remove zero uncertainties
    good_data = data_input[np.where(data_input[:, 2] != 0)]

    # remove Nans and infinities
    for index in range(len(data_input[0])):
        good_data = good_data[np.where(np.isfinite(good_data[:, index]))]

    # remove obvious outliers
    error_indicies = np.array([], dtype='int')

    for index in range(len(good_data[:, 0]) - 2):

        upper_difference = np.abs(good_data[index, 1] - good_data[index + 1, 1])
        lower_difference = np.abs(good_data[index + 2, 1] - good_data[index + 1, 1])
        allowed_gap = good_data[index + 1, 2] * 10

        if upper_difference > allowed_gap and lower_difference > allowed_gap:
            error_indicies = np.append(error_indicies, index + 1)

    good_data = np.delete(good_data, error_indicies, axis=0)
    return good_data


def remove_outliers(data_input, model, model_parameters):
    """
    Removes outliers from the data, which are more that 5 standard deviations
    (error value) away from the model.

    Parameters
    ----------
    data_inut : N-dimentional array of floats
    model : function
    model perameters : array of floats

    Returns
    -------
    correted_data : N-dimentional array of floats

    """

    value_difference = np.abs(data_input[:, 1] - model(model_parameters,
                                                       data_input[:, 0] * pc.year))
    correted_data = data_input[np.where(data_input[:, 2] * 5 > value_difference)]

    return correted_data


def star_velocity_data(wavelengths, wavelength_uncertainties, rest_wavelenght,
                       inclination_angle=np.pi/2):
    """
    Converts Doppler shift wavelengts data into velocoties data.

    Parameters
    ----------
    wavelengths : float array
    wavelengths_uncertainties : float array
    rest_wavelenght : float
    inclination_angle : float
        Angle between the line of sight of the observer and the normal to the
        planet's orbital plane.

    Returns
    -------
    velocoties: float array
    velocity_uncertainties: float array

    """
    velocities = (pc.speed_of_light / np.sin(inclination_angle)
                                    * (1 - wavelengths / rest_wavelenght))
    velocity_uncertainties = np.abs(pc.speed_of_light
                                    / (rest_wavelenght
                                       * np.sin(inclination_angle))
                                    * wavelength_uncertainties)

    return velocities, velocity_uncertainties


def star_velocity_model(velocity_magnitude_and_angular_velocity, time,
                        inclination_angle=np.pi/2, initial_phase=0):
    """
    The expected model for the velocty variation of a star.
    v_s(t) = v_0 sin(ωt + φ)

    Parameters
    ----------
    time : float
    velocity_magnitude_and_angular_velocity : float array
        has two values within it: magnitude of the star's velocity and its
        angular velocity
    inclination_angle : float
        Angle from the line of sight of the observer and the normat to the
        planet's orbital plane.
    initial_phase : float
        The default is 0

    Returns
    -------
    star_velocity : float

    """
    velocity_magnitude = velocity_magnitude_and_angular_velocity[0]
    angular_momentum = velocity_magnitude_and_angular_velocity[1]
    star_velocity = (velocity_magnitude * (np.sin(angular_momentum * time
                                                  + initial_phase))
                     / np.sin(inclination_angle))

    return star_velocity


def chi_squared(prediction, data, uncertainty):
    """
    Calculates the value of chi squared for the fit.

    Parameters
    ----------
    prediction : array of floats
    data : array of floats
    uncertainty : array of floats

    Returns
    -------
    Chi squared : float

    """

    return np.sum(((prediction - data) / uncertainty)**2)


def planet_mass(star_mass, star_velocity, planet_velocity_value):
    """
    Calculates the mass of the planet given parameters, specified below.
    m_planet = M_sstar * v_star / v_planet

    Parameters
    ----------
    star_mass : float
    star_velocity : float
    planet_velocity : float

    Returns
    -------
    planet_mass_value : float

    """
    planet_mass_value = star_mass * star_velocity / planet_velocity_value
    return planet_mass_value


def planet_mass_uncertainty(star_velocity, star_velocity_uncertainty,
                            planet_velocity_value, planet_velocity_uncertainty,
                            planet_mass_value):
    """
    Calculates the uncertainty on the mass of the planet defined above.

    Parameters
    ----------
    star_velocity : float
    star_velocity_uncertainty : float
    planet_velocity_value : float
    planet_velocity_uncertainty : float
    planet_mass_value : float

    Returns
    -------
    planet_mass_uncertainty : float

    """
    planet_mass_uncertainty = np.sqrt((star_velocity_uncertainty
                                       / star_velocity)**2
                                      + (planet_velocity_uncertainty
                                         / planet_velocity_value)**2) * planet_mass_value
    return planet_mass_uncertainty


def planet_velocity(star_mass, orbital_distance_value):
    """
    Calculates orbital velocity of the planet given parameters, specified below.
    v_planet = (G * M_star / r)^1/2

    Parameters
    ----------
    star_mass : float
    orbital_distance_value : float

    Returns
    -------
    planet_velocity_value : float

    """
    planet_velocity_value = np.sqrt(pc.gravitational_constant * star_mass
                                    / orbital_distance_value)
    return planet_velocity_value


def planet_velocity_uncertainty(star_mass, orbital_distance_value,
                                orbital_distane_uncertainty):
    """
    Calculates the uncertainty on the velocity of the platet defined above.

    Parameters
    ----------
    star_mass : float
    orbital_distance_value : float
    orbital_distane_uncertainty : float

    Returns
    -------
    planet_velocity_uncertainty : float

    """
    planet_velocity_uncertainty = (0.5 * np.sqrt(pc.gravitational_constant
                                                 * star_mass)
                                   * orbital_distance_value**(-3 / 2)
                                   * orbital_distane_uncertainty)
    return planet_velocity_uncertainty


def orbital_distance(star_mass, angular_velocity):
    """
    Calculates orbital distance r from the planet to the star
    given parameters, specified below.

    r = (G * M_star * P^2 / 4 * pi^2)^1/3
    P = 2pi / ω

    Parameters
    ----------
    star_mass : float
    angular_velocity : float

    Returns
    -------
    orbital_distance_value : foat

    """
    orbital_distance_value = np.cbrt(pc.gravitational_constant * star_mass
                                     * (2 * np.pi / angular_velocity)**2
                                     / (4 * np.pi**2))
    return orbital_distance_value


def orbital_distance_uncertainty(star_mass, angular_velocity,
                                 angular_velocity_uncertainty):
    """
    Calculates the uncertainty on the orbital distance defined above.

    Parameters
    ----------
    star_mass : float
    angular_velocity : float
    angular_velocity_uncertainty : float

    Returns
    -------
    orbital_distance_uncertainty : float

    """
    orbital_distance_uncertainty = ((2 / 3) * np.cbrt(pc.gravitational_constant
                                                      * star_mass)
                                    * angular_velocity**(-5 / 3)
                                    * angular_velocity_uncertainty)
    return orbital_distance_uncertainty


def reduced_chi_squared(chi_squared_value, data, parameter_number):
    """
    Calculates reduced chi squared from a full chi squared of the data fit.
    reduced_chi^2 = chi^2 / Ndof
    Ndof is the number od degrees of freeedom.

    Parameters
    ----------
    chi_squared : float
    data_set : array of floats
    parameter_number : float

    Returns
    -------
    reduced_chi_squared_value : float

    """

    reduced_chi_squared_value = chi_squared_value / (len(data) - parameter_number)
    return reduced_chi_squared_value


def model_plot(x_values, y_values, y_errors, model_function, model_parameters):
    """
    Makes a plot of data with its uncertainties and model to compare with.

    Parameters
    ----------
    x_values : float array
        Should represent time in years.
    y-values : float array
        Should represent velocity on kilometers per second.
    y_errors : float array
    model_function : function
    model_parameters : float array

    Returns
    -------
    0

    """

    figure = plt.figure()
    axes = figure.add_subplot(111)
    axes.scatter(x_values, y_values, y_errors, label='data')
    axes.errorbar(x_values, y_values, y_errors, linestyle='')
    axes.plot(x_values, model_function(model_parameters, x_values * pc.year),
              label='model')
    axes.set_title('Velocity of the Star Obtained Using Redshift')
    axes.set_xlabel('time / years')
    axes.set_ylabel('velocity / kms$^{-1}$')
    axes.legend()
    plt.savefig('velocity_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

    return 0


def mesh_arrays(x_array, y_array):
    """
    Returns two meshed arrays of size len(x_array)
    by len(y_array)

    Parameters
    ----------
    x_array : float array
    y_array : float array

    Returns
    -------
    x_array_mesh : 2D float array
    y_array_mesh : 2D float array

    """
    x_array_mesh = np.empty((0, len(x_array)))

    for dummy_element in y_array:
        x_array_mesh = np.vstack((x_array_mesh, x_array))

    y_array_mesh = np.empty((0, len(y_array)))

    for dummy_element in x_array:
        y_array_mesh = np.vstack((y_array_mesh, y_array))

    y_array_mesh = np.transpose(y_array_mesh)

    return x_array_mesh, y_array_mesh


def contour_plot(data, model, model_parameters):
    """
    Makes a contour plot of chi squared values around a minimum with 2
    parameters varied.

    Parameters
    ----------
    data : float array
        must have shape (n, 3)
    model : function
    model_parameters : float array
        2 optimal parameters for the fit

    Returns
    -------
    None.

    """

    x_values = np.linspace(-25, 75, 100)
    y_values = np.linspace(1.5e-8, 4e-8, 100)

    chi_squares = np.empty([100, 100])

    for index_x, value_x in enumerate(x_values):
        for index_y, value_y in enumerate(y_values):
            chi_squares[index_y, index_x] = chi_squared(model([value_x, value_y],
                                                              data[:, 0] * pc.year),
                                                        data[:, 1],
                                                        data[:, 2])

    x_values_mesh, y_values_mesh = mesh_arrays(x_values, y_values)
    figure = plt.figure(figsize=(4, 4))
    axis = figure.add_subplot(111)
    axis.contour(x_values_mesh, y_values_mesh, chi_squares, 50)
    axis.scatter(model_parameters[0], model_parameters[1],
                 color='k', label='Minimum')
    axis.legend()
    axis.set_title(r'$\chi^2$ contour plot', fontsize=14)
    axis.set_xlabel('maximum velocity / kms$^{-1}$')
    axis.set_ylabel('angular velocity / rads$^{-1}$')
    # axis.clabel(CONTOUR_PLOT)
    plt.savefig('contour_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

    return 0


def uncertainty_contour_plot(data, model, model_parameters, chi_squared_min):
    """
    Makes a contour plot of chi squared values around a minimum with 2
    parameters varied. Calculates the std for the fit parameters.

    Parameters
    ----------
    data : float array
        must have shape (n, 3)
    model : function
    model_parameters : float array
        2 optimal parameters for the fit

    Returns
    -------
    velocity_uncertainty : float
    angular_velocity_uncertainty : float

    """
    # creating a plot
    x_values = np.linspace(27, 33, 100)
    y_values = np.linspace(2.35e-8, 2.6e-8, 100)

    chi_squares = np.empty([100, 100])

    for index_x, x_value in enumerate(x_values):
        for index_y, y_value in enumerate(y_values):
            chi_squares[index_y, index_x] = chi_squared(model([x_value, y_value],
                                                              data[:, 0] * pc.year),
                                                        data[:, 1],
                                                        data[:, 2])

    x_values_mesh, y_values_mesh = mesh_arrays(x_values, y_values)
    uncertainty_figure = plt.figure()
    chi_squared_levels = [chi_squared_min + 1, chi_squared_min + 2.30,
                          chi_squared_min + 5.99, chi_squared_min + 9.21]

    uncertainty_axis = uncertainty_figure.add_subplot(111)
    uncertainty_plot = uncertainty_axis.contour(x_values_mesh, y_values_mesh,
                                                chi_squares,
                                                levels=chi_squared_levels)
    uncertainty_axis.scatter(model_parameters[0], model_parameters[1],
                             color='k')
    uncertainty_axis.set_title(r'       $\chi^2$ contour uncertainties plot',
                               fontsize=14)
    uncertainty_axis.set_xlabel('maximum velocity / kms$^{-1}$')
    uncertainty_axis.set_ylabel('angular velocity / rads$^{-1}$')
    labels = [r'$\chi^2_{{\mathrm{{min.}}}}+1.00$',
              r'$\chi^2_{{\mathrm{{min.}}}}+2.30$',
              r'$\chi^2_{{\mathrm{{min.}}}}+5.99$',
              r'$\chi^2_{{\mathrm{{min.}}}}+9.21$',
              'Minimum']

    uncertainty_axis.clabel(uncertainty_plot)
    box = uncertainty_axis.get_position()
    uncertainty_axis.set_position([box.x0, box.y0, box.width * 0.7,
                                   box.height])

    # Add custom plot labels
    for index, label in enumerate(labels):
        uncertainty_axis.collections[index].set_label(label)
    uncertainty_axis.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                            fontsize=12)
    plt.savefig('uncertainty_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

    # finding the uncertainties
    data_chi = uncertainty_plot.allsegs[0][0]
    velocity_uncertainty = (np.max(data_chi[:, 0]) - np.min(data_chi[:, 0])) / 2
    angular_velocity_uncertainty = (np.max(data_chi[:, 1]) - np.min(data_chi[:, 1])) / 2

    return velocity_uncertainty, angular_velocity_uncertainty


def reduced_chi_squared_check(chi_squared_reduced):
    """
    Checks if the reduced chi squared lies in the allowed range of
    0.5 inclusive to 1.5 not inclusive. Prints a warning message if the
    reduced chi squared is not suitable for a good fit.

    Parameters
    ----------
    reduced_chi_squared : float

    Returns
    -------
    0

    """
    if chi_squared_reduced <= 0.5 or chi_squared_reduced > 1.5:
        print('Warning! Reduced chi squared is out of allowed range. Fit might '
              'not represent the data well.')
    return 0


def find_planet_parameters(star_mass, optimal_parameters):
    """
    Calculates planet's orbital distance and mass given the mass of the star it
    orbits and the parameters for the velocity fit.

    Parameters
    ----------
    star_mass : float
    optimal_parameters : float

    Returns
    -------
    planets_orbital_distance : float
    planets_mass : float

    """

    planets_orbital_distance = orbital_distance(star_mass, optimal_parameters[1])
    planets_mass = planet_mass(star_mass, optimal_parameters[0],
                               planet_velocity(star_mass, planets_orbital_distance))
    return planets_orbital_distance, planets_mass


def main(file_names, rest_wavelength, star_mass, jovian_mass):
    """
    Main code. Prints out the required planet and fit parameters, and creates
    some plots representing the data and the fit.

    Parameters
    ----------
    file_names : float
    rest_wavelength : float
    star_mass : float
    jovian_mass : float

    Returns
    -------
    None

    """

    raw_data, file_found = read_and_combine(file_names)

    if file_found:
        # validating data and removing outliers
        initial_data = validate(raw_data)
        velocities, velocity_uncertainties = star_velocity_data(initial_data[:, 1],
                                                                initial_data[:, 2],
                                                                rest_wavelength)
        optimal_parameters = fmin(lambda x: chi_squared(star_velocity_model(x,
                                                                            (initial_data[:, 0]
                                                                             * pc.year)),
                                                        velocities,
                                                        velocity_uncertainties),
                                  [50, 3e-8], disp=0)
        initial_data_velocity = np.transpose(np.vstack((initial_data[:, 0], velocities,
                                                        velocity_uncertainties)))
        correted_data = remove_outliers(initial_data_velocity,
                                        star_velocity_model,
                                        optimal_parameters)
        velocities, velocity_uncertainties = correted_data[:, 1], correted_data[:, 2]

        # finding the fit parameters
        optimal_parameters = fmin(lambda x: chi_squared(star_velocity_model(x,
                                                                            (correted_data[:, 0]
                                                                             * pc.year)),
                                                        velocities,
                                                        velocity_uncertainties),
                                  [50, 3e-8])
        chi_squared_min = chi_squared(star_velocity_model(optimal_parameters,
                                                          correted_data[:, 0] * pc.year),
                                      velocities,
                                      velocity_uncertainties)

        reduced_chi_squared_value = reduced_chi_squared(chi_squared_min, velocities, 2)
        reduced_chi_squared_check(reduced_chi_squared_value)

        # creating plots

        model_plot(correted_data[:, 0], velocities, velocity_uncertainties,
                   star_velocity_model, optimal_parameters)
        contour_plot(correted_data, star_velocity_model, optimal_parameters)
        velocity_uncertainty, angular_velocity_uncertainty = uncertainty_contour_plot(correted_data,
                                                                                      star_velocity_model,
                                                                                      optimal_parameters,
                                                                                      chi_squared_min)

        # finding planet parameters

        orbit_radius, planets_mass = find_planet_parameters(star_mass,
                                                            optimal_parameters)
        radius_uncertainty = orbital_distance_uncertainty(star_mass,
                                                          optimal_parameters[1],
                                                          angular_velocity_uncertainty)
        mass_uncertainty = planet_mass_uncertainty(optimal_parameters[0],
                                                   velocity_uncertainty,
                                                   planet_velocity(star_mass,
                                                                   orbit_radius),
                                                   planet_velocity_uncertainty(star_mass,
                                                                               orbit_radius,
                                                                               radius_uncertainty),
                                                   planets_mass)

        # print the results
        print('')
        print('Fit parameters:')
        print('max velocity of the star v_0 = {:2.2f} +/- {:2.2f} ms^-1'
              .format(optimal_parameters[0], velocity_uncertainty))
        print('angular velocity of the star ω = {:1.3e} +/- {:1.1e} rads^-1'
              .format(optimal_parameters[1], angular_velocity_uncertainty))
        print('reduced chi^2 = {:0.3f}'.format(reduced_chi_squared_value))
        print('')
        print('Planet parameters:')
        print('mass m_p = {:1.3f} +/- {:1.3f} Jovian masses'
              .format(planets_mass / jovian_mass, mass_uncertainty / jovian_mass))
        print('orbital distane r = {:1.3f} +/- {:1.3f} AU'
              .format(orbit_radius / pc.astronomical_unit,
                      radius_uncertainty / pc.astronomical_unit))

main(FILE_NAMES, REST_WAVELENGHT, STAR_MASS, JOVIAN_MASS)
