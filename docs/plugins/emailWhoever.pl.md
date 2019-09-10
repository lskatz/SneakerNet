# SYNOPSIS

`emailWhoever.pl`

Emails all SneakerNet results

# Software requirements

* sendmail (already installed on most systems)

# Algorithm

Attaches all files under `SneakerNet/forEmail` and emails
them to the list under `config/emails.conf`. If any emails
are found in `snok.txt` or `SampleSheet.csv`, then they 
will also receive the email.

# Outputs

email with attachments

