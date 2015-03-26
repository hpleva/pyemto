Installation
============

Python in an interpreted language, therefore nothing needs to be compiled.
Follow these steps to install pyEMTO:

1. Copy pyEMTO to a folder of your choice.
2. Open ~/.bashrc in a text editor and add the following lines:

   .. sourcecode:: bash

      PYTHONPATH="${PYTHONPATH}:/path/to/pyEMTO"
      export PYTHONPATH

3. Save and close.
4. **Optional**: The new path can be activated without restarting your terminal session by executing:

   .. sourcecode:: bash

      source ~/.bashrc

5. **Optional**: Check the installation process has succeeded by typing:

   .. sourcecode:: bash

      python
      
   .. sourcecode:: bash

      >>>import sys

   .. sourcecode:: bash

      >>>print(sys.path)

   If '/path/to/pyEMTO' can be found within the printout, PYTHONPATH has been succesfully set. Finally, type:

   .. sourcecode:: bash

      >>>import pyemto

   followed by

   .. sourcecode:: bash

      >>>help(pyemto)

   to verify that pyEMTO has been installed succesfully.

6. **Optional**: Try running some of the scripts from the :ref:`examples <examples-page>` page to verify the integrity of your pyEMTO installation.
