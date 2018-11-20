
Implementation of TITE-PK 
=========================
This is the accompanying code to paper "Phase I dose-escalation trials with more than one dosing regimen"
by BK Günhan, S Weber, A Seroutou, and T Friede.
See the Preprint version: 

TITE-PK 
=======
TITE-PK (a time-to-event pharmacokinetic model) is a method for designing and analyzing 
phase I dose-escalation trials with dose-regimen changes, eg a daily and weekly dosing.
Hence, it is a model-based approach for dose-escalation trials like the Bayesian Logistic
Regression Model, but TITE-PK is able to integrate varying dosing regimens. For details, 
see the preprint.

Computations
============
TITE-PK uses **Stan** via `Rstan`, hence it should be installed. 
See [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
for how to install it.

`TITEPK_run.R` is the main R-script which produce the Figure 2 of the main text
for TITE-PK methods.

`tite_pk_one_param.stan` and `tite_pk_two_params.stan` are main **Stan** scripts
which are used in `TITEPK_run.R` to implement the TITE-PK models.


Licence
=======

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.




