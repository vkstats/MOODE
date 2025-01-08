## Resubmission

This is a resubmission. 

* Search.R: replaced the call to `print` with a call to `cli::cli_alert_info`
and renamed the `info.update` argument to `verbose`.

* print.settings.R: added a `\value` tag specifying that no values are returned.
A description has been added to the help file, indicating the parameters which
are summarised.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
