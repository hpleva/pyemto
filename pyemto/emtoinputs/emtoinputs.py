# -*- coding: utf-8 -*-
"""

Created on Wed Dec  3 14:25:06 2014

@author: Matti Ropo
@author: Henrik Levämäki
"""


class Emtoinputs:

    """Class which is used to communicate with the Kgrn and Kfcd classes.

    :returns: None
    :rtype: None
    """

    def __init__(self):

        # Import necessary packages
        from pyemto.emtoinputs.kgrn import Kgrn
        from pyemto.emtoinputs.kfcd import Kfcd
        from pyemto.emtoinputs.batch import Batch

        self.kgrn = Kgrn()
        self.kfcd = Kfcd()
        self.batch = Batch()

        return

    def set_values(self, **kwargs):
        """Passes various input parameters down to the Kgrn and Kfcd classes

        :param **kwargs: Keyword arguments
        :type **kwargs: dict
        :returns: None
        :rtype: None
        """

        for key in kwargs:
            attr_found = False
            if hasattr(self.kgrn, key):
                self.kgrn.set_values(key, kwargs[key])
                attr_found = True
            if hasattr(self.kfcd, key):
                self.kfcd.set_values(key, kwargs[key])
                attr_found = True
            if hasattr(self.batch, key):
                self.batch.set_values(key, kwargs[key])
                attr_found = True
            # if attr_found == False:
                # print('WARNING: Neither Kgrn(), Kfcd() nor Batch_emto() classes' +
                      # ' have the attribute \'{0}\''.format(key))
        return
