gravica
=======

**gravica** is a General Relativity computation library built on
`Symbolica <https://symbolica.io>`_.

It provides a pipeline of tensor classes that lazily compute differential-geometry
objects from a metric tensor:

.. code-block:: text

   MetricTensor → ChristoffelSymbols → RiemannTensor → RicciTensor → EinsteinTensor / WeylTensor
                         ↓                   ↓              ↓              ↓
                  GeodesicEquations   KretschnerScalar  SchoutenTensor  StressEnergyTensor

.. toctree::
   :maxdepth: 2
   :caption: Contents

   api/index
