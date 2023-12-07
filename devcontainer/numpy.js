{{! Numpy Docstring Template }}
{{summaryPlaceholder}}

{{extendedSummaryPlaceholder}}
{{#parametersExist}}

INPUT:

{{#args}}

- ``{{var}}`` -- {{typePlaceholder}}  {{descriptionPlaceholder}}
{{/args}}
{{#kwargs}}

- ``{{var}}`` : {{typePlaceholder}}, optional  {{descriptionPlaceholder}}, by default {{&default}}
{{/kwargs}}
{{/parametersExist}}

{{#returnsExist}}

OUTPUT:

{{#returns}}
{{typePlaceholder}}  {{descriptionPlaceholder}}
{{/returns}}
{{/returnsExist}}


{{#yieldsExist}}
YIELD_OUTPUT:

``{{#yields}}``
{{typePlaceholder}}  {{descriptionPlaceholder}}
{{/yields}}
{{/yieldsExist}}
{{#exceptionsExist}}

EXCEPTIONS:
----------
{{#exceptions}}
{{type}}  {{descriptionPlaceholder}}
{{/exceptions}}
{{/exceptionsExist}}

EXAMPLES::

