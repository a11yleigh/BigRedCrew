library(gert)
library(credentials)
gert::git_config_global_set("user.name", "rutved118")
gert::git_config_global_set("user.email", "rutved118@gmail.com")   # use the email on your GitHub account

credentials::git_credential_forget("https://github.com")
credentials::set_github_pat()   # paste the NEW token
# prompts a popup asking you to enter your GitHub Personal Access Token

gert::git_pull()                          # pull the most recent changes from GitHub
gert::git_add(dir(all.files = TRUE))      # stage every new or edited file
gert::git_commit_all("my first commit")   # save your record of the edits (a commit)
gert::git_push()                          # push your commit up to GitHub

