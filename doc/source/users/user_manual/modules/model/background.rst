Background models
=================

Background models differ from sky models in that those models are not convolved
with the instrument response function, and directly provide the number of
predicted background counts as function of measured event properties.
Background models need to be implemented in the instrument specific models.
All background models derive from the abstract :doxy:`GModelData` base
class.



