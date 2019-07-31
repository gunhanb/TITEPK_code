Implementation of TITE-PK 
=========================
This is the accompanying code to paper "A Bayesian time-to-event pharmacokinetic model for phase I
dose-escalation trials with multiple schedules"
by BK Günhan, S Weber, and T Friede.
We consider an artificial phase I trial with multiple schedules to illustrate our proposed model TITE-PK.

TITE-PK 
=======
TITE-PK (a time-to-event pharmacokinetic model) is a method for designing and analyzing 
phase I dose-escalation trials with multiple schedules, eg a daily and weekly dosing.
Hence, TITE-PK finds an acceptable dose and schedule combination. 

Computations
============
TITE-PK uses **Stan** via `Rstan`, hence it should be installed. 
See [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
for how to install it.

`TITEPK_run.R` is the main R-script which produce the Figure 2 of the main text
for TITE-PK methods.

`tite_pk.stan` is the main **Stan** script
which is used in `TITEPK_run.R` to implement the TITE-PK model.

Session info
============
```{r, results='markup'} 
sessionInfo()
```



Licence
=======

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  The code was written by [Burak Kürsad
  Günhan](http://bkguenhan.rbind.io) and Sebastian Weber. To report any
  issues or bugs or to suggest enhancements to the package, please go
  [here](https://github.com/gunhanb/TITEPK_code/issues).




