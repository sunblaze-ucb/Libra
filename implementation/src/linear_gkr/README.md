# Interface for Prover & Temp Interface for Verifier

## File Format

In the first line, input a integer `d` represent the number of layers.
Input layer has layer number `0` output layer has layer number `d-1`

For next `d - 1` lines, each line specify a layer.

For `i`-th line, it specify the layer `i`, the first number is `n`, specify the number of gates in this layer. `n` must be a multiple of `2`.
The rest of this line contains `4n` integers, represent `n` gates. For each gate, we use for integers to describe: `ty g u v`, indicates the type of the gate, and the connection of the gate, `g` is the gate number, `u` is the left input of the gate, `v` is the right input of the gate.

We have `3` different types of gates for now. They are addition, multiplication, dummy. Dummy gate takes no input, and never serves as a input of other gate.

`ty=0` is addition gate, `ty=1` is multiplication gate, `ty=2` is dummy gate.

## Example
``
3 \\two layers
2 0 0 0 1 1 1 2 3 \\first gate is addition, and second is a multiplication
1 1 0 0 1 \\this is the output layer, it's a multiplication gate
``