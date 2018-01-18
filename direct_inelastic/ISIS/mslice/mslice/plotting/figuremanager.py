# This is the class responsible for handling which figures are current (will receive plotting commands) and which are
# kept (will not be modified until it is made current once again).
# The FigureManager in its responsibilities highly resembles the class maptplotlib._pylab_helpers.Gcf, However it
# Adds the functionality of having multiple 'categories' which category having it own current window
# This is achieved through use of the supplied decorator 'activate_category' This decorator accepts one parameter  (a
# string ) specifying which category a function belongs two. For instance to apply the the the decorator
# activate_category('<category>') to the the function pyplot.pcolor would signal that the function pcolor would only apply
# to  plots of the category '<category>'. All of this is done through manipulating the the return value of `pyplot.gcf`
# (If you are not familiar with with how matplotlib.pyplot functions work now would be a good time to go and familiarize
# yourself. Most of the functions are autogenerated and are almost indentical, looking the gca, gcf, plot and imshow
# functions should be enough to get yourself up to speed.)
#
# gcf in pyplot should return to FigureManager.get_active_figure()
#
# If gcf is called from inside a categorized (decorated) function then it will return the current figure for that
# functions category. If there is no current figure for that category then it will create a new figure and return it.
#
# If a new figure is created in a call of gcf that is categorized then the new figure will automatically be assigned
# to the category. And it will be set as the current figure for that category.
#
# If gcf is called from inside an uncategorized function then it should return the `active figure` which is defined as
# the last plot window to receive any plotting command (regardless of category). This makes sense for functions which
# can apply to any plots such as `xlabel`
#
# if a new figure is created by an uncategorized function then will be 'uncategorized'. If the current 'active figure'
# is an uncategorized figure and categorized function is called then it should be returned and then that figure should
# be added to the category of the command
#
# Currently there are only two categories ('1d' and '2d') Hard-coded into the manager however it would be simple to add
#  more by simply adding more entries into the _figures_by_category and _category_current_figures
#  dictionaries (Or maybe adding an add_category function to facilitate this)
#
#
from __future__ import (absolute_import, division, print_function)
from functools import wraps

from mslice.plotting import get_figure_class


class NoFigure(object):
    def __init__(self):
        pass

    def __eq__(self, other):
        return isinstance(other, NoFigure)

    def __getattribute__(self, item):
        if item in ('__init__', '__eq__'):
            return super(NoFigure, self).__getattribute__(item)
        raise RuntimeError('Attempt to manipulate non-existent figure')


PlotFigureManager = get_figure_class()
NO_FIGURE = NoFigure()


class FigureManager(object):
    """This is singleton static class to manage the current _figures

    Do not instantiate this class"""
    # if there is a current figure it should be both current and active
    _active_category = None
    _category_current_figures = {"1d": NO_FIGURE, "2d": NO_FIGURE}  # Current _figures recieve decorated commands
    _figures_by_category = {"1d": [], "2d": []}
    _unclassified_figures = []
    _active_figure = NO_FIGURE  # Will receive all commands that have a matching decorator or are undecorated
    _figures = {}

    def __init__(self, *args, **kwargs):
        raise Exception("This is a static class singleton. Do not Instantiate it")

    @staticmethod
    def _new_figure(fig_num=None):
        if fig_num is None:
            fig_num = 1
            while any([fig_num == existing_fig_num for existing_fig_num in FigureManager._figures.keys()]):
                fig_num += 1
        new_fig = PlotFigureManager(fig_num, FigureManager)
        FigureManager._figures[fig_num] = new_fig
        FigureManager._active_figure = fig_num
        FigureManager._unclassified_figures.append(fig_num)
        return FigureManager._figures[fig_num], fig_num

    @staticmethod
    def get_figure_number(fig_num=None, create_if_not_found=True):
        FigureManager._active_figure = fig_num
        try:
            # make the figure the current figure of its category
            category = FigureManager.get_category(fig_num)
            FigureManager._category_current_figures[category] = fig_num
        except KeyError:
            # the figure is still uncategorised, Do nothing
            pass

        figure = FigureManager._figures.get(fig_num, None)
        if figure or not create_if_not_found:
            return figure
        else:
            # return the figure discard the number
            return FigureManager._new_figure(fig_num)[0]


    @staticmethod
    def get_active_figure():
        if FigureManager._active_category:
            if FigureManager._active_figure in FigureManager._unclassified_figures:
                FigureManager.assign_figure_to_category(FigureManager._active_figure, FigureManager._active_category,
                                                        make_current=True)


            elif FigureManager._category_current_figures[FigureManager._active_category] == NO_FIGURE:
                _, num = FigureManager._new_figure()
                FigureManager.assign_figure_to_category(num, FigureManager._active_category, make_current=True)
                FigureManager._active_figure = num

            else:
                FigureManager._active_figure = FigureManager._category_current_figures[FigureManager._active_category]

        else:
            if FigureManager._active_figure == NO_FIGURE:
                fig, num = FigureManager._new_figure()
                FigureManager._active_figure = num

        return FigureManager._figures[FigureManager._active_figure]


    @staticmethod
    def _activate_category(category):
        """Sets the active category to the supplied argument, do not call this function directly. Instead use supplied
        Decorator below 'activate_category' """
        FigureManager._active_category = category

    @staticmethod
    def _deactivate_category():
        """ Unsets the active category. do not call this function directly. Instead use supplied
        Decorator below 'activate_category' """
        FigureManager._active_category = None

    @staticmethod
    def assign_figure_to_category(fig_num, category, make_current=False):
        if fig_num not in FigureManager._figures:
            raise ValueError("Figure does not exist")

        if fig_num in FigureManager._unclassified_figures:
            FigureManager._unclassified_figures.remove(fig_num)

        for a_category in FigureManager._figures_by_category:
            if fig_num in FigureManager._figures_by_category[a_category]:
                FigureManager._figures_by_category[a_category].remove(fig_num)
            if FigureManager._category_current_figures == fig_num:
                FigureManager._category_current_figures = NO_FIGURE

        FigureManager._figures_by_category[category].append(fig_num)
        if make_current:
            FigureManager._category_current_figures[category] = fig_num
        FigureManager.broadcast()

    @staticmethod
    def figure_closed(figure_number):
        """Figure is closed, remove all references to it from all internal list

        If it was the category current or global active figure then set that to NO_FIGURE"""
        if FigureManager._active_figure == figure_number:
            FigureManager._active_figure = NO_FIGURE
        for a_category in FigureManager._figures_by_category:
            if figure_number in FigureManager._figures_by_category[a_category]:
                FigureManager._figures_by_category[a_category].remove(figure_number)

            if FigureManager._category_current_figures[a_category] == figure_number:
                FigureManager._category_current_figures[a_category] = NO_FIGURE
        try:
            del FigureManager._figures[figure_number]
        except KeyError:
            raise KeyError('The key "%s" does not exist. The figure cannot be closed' % figure_number)

    @staticmethod
    def get_category(figure_number):
        """Return the category of the figure"""
        for category,fig_list in list(FigureManager._figures_by_category.items()):
            if figure_number in fig_list:
                figure_category = category
                break
        else:
            raise KeyError("Figure no. %i was not found in any category "%figure_number if figure_number else 0)
            # in-line if handles the case figure_number is None
        return figure_category

    @staticmethod
    def set_figure_as_kept(figure_number):
        # kept figures are just lying around, not really managed much, until they report in as current again
        if FigureManager._active_figure == figure_number:
            FigureManager._active_figure = NO_FIGURE
        try:
            figure_category = FigureManager.get_category(figure_number)
        except KeyError:
            figure_category = None

        if figure_category:
            if FigureManager._category_current_figures[figure_category] == figure_number:
                FigureManager._category_current_figures[figure_category] = NO_FIGURE

        FigureManager.broadcast(figure_category)

    @staticmethod
    def set_figure_as_current(figure_number):
        try:
            figure_category = FigureManager.get_category(figure_number)
        except KeyError:
            figure_category = None
        if figure_category:
            FigureManager._category_current_figures[figure_category] = figure_number
        FigureManager._active_figure = figure_number
        FigureManager.broadcast(figure_category)


    @staticmethod
    def broadcast(category=None):
        """This method will broadcast to all figures in 'category' to update the displayed kept/current status"""
        if category is None:
            broadcast_list = FigureManager._figures_by_category
        else:
            broadcast_list = [category]

        for category in broadcast_list:
            for figure_number in FigureManager._figures_by_category[category]:

                if FigureManager._category_current_figures[category] == figure_number:
                    FigureManager._figures[figure_number].set_as_current()

                else:
                    FigureManager._figures[figure_number].set_as_kept()

        for figure in FigureManager._unclassified_figures:
            if figure == FigureManager._active_figure:
                FigureManager._figures[figure].set_as_current()
            else:
                FigureManager._figures[figure].set_as_kept()

    @staticmethod
    def all_figure_numbers():
        """An iterator over all figure numbers"""
        return list(FigureManager._figures.keys())

    @staticmethod
    def all_figures_numbers_in_category(category):
        """Return an iterator over all _figures numbers in a category"""
        return iter(FigureManager._figures_by_category[category])

    @staticmethod
    def unclassified_figures():
        """Return an iterator over all unclassified _figures"""
        return iter(FigureManager._unclassified_figures)

    @staticmethod
    def reset():
        """Reset all class variables to initial state. This Function exists for testing purposes """

        FigureManager._active_category = None
        FigureManager._category_current_figures = {"1d": NO_FIGURE, "2d": NO_FIGURE}  # Current _figures are overplotted
        FigureManager._figures_by_category = {"1d": [], "2d": []}
        FigureManager._unclassified_figures = []
        FigureManager._active_figure = NO_FIGURE
        FigureManager._figures = {}

    @staticmethod
    def all_figures():
        """Return an iterator over all figures"""
        return list(FigureManager._figures.values())

    @staticmethod
    def number_of_figure(figure):
        for key,value in list(FigureManager._figures.items()):
            if value == figure:
                return key
        raise ValueError('Figure %s was not recognised'%figure)


def activate_category(category):

    def real_activate_function_decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            FigureManager._activate_category(category)
            return_value = function(*args, **kwargs)
            FigureManager._deactivate_category()
            return return_value
        return wrapper

    return real_activate_function_decorator
