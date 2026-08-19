****************
[backend-metal]
****************

Parameterises the Metal backend with

#. ``gimmik-max-nnz`` --- cutoff for GiMMiK in terms of the number of
   non-zero entires in a constant matrix, overriding the suitability
   criteria GiMMiK applies for the GPU family in question:

    *int*

#. ``gimmik-nkerns`` --- number of kernel algorithms to try when
   benchmarking, defaults to 18:

    *int*

#. ``gimmik-nbench`` --- number of benchmarking runs for each
   kernel, defaults to 40:

     *int*

Example:

.. code-block:: ini

    [backend-metal]
    gimmik-max-nnz = 512
