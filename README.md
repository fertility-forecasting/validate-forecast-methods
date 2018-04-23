# validate-forecast-methods

`validate-forecast-methods` provides the R source code used in the project: 

_Bohk-Ewald, Christina, Peng Li, and Mikko Myrskylä (2017). Assessing the accuracy of cohort fertility forecasts. Presented in session: Statistical methods in demography at the PAA 2017 Annual Meeting, Chicago, IL, USA, April 27-April 29, 2017_

## How to run

The code needs to be executed in eight steps, as defined by `step-*.R` files in the root directory. You can run these scripts separately or all at once with `step-0-all-steps-at-once.R` 

### Prerequisites

- The `.txt` files in `./input-data/` contain empty lists due to copyrights. You will need to put in fertility data yourself. See comments for step 2 in [`step-0-all-steps-at-once.R`].

- Make sure you have configured the correct paths (in each of the nine `step-*.R` files).

### Execution

You can run all steps all at once `step-0-all-steps-at-once.R` (recommended) or you can run the `step-*.R` scripts in the root directory in the prescribed order.

## How to cite

If you use this code for academic research, please cite the above conference talk.

## How to contribute

Please note that this source code is an academic project. We welcome any issues and pull requests.

## License

The source code of `validate-forecast-methods` is published under the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). 