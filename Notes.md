# Notes

### Update R version (3.4 -> 3.5) in Ubuntu 16.04

```
sudo su
echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
apt-get update
apt-get install r-base
apt-get install r-base-dev
```
Remember to check the Ubuntu version name using `sudo lsb_release -a`

See [CRAN](https://cran.r-project.org/bin/linux/ubuntu/) also.

### Locale setting
```
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"
sudo dpkg-reconfigure locales
```

### Download AWS s3 bucket objects

- Install AWS CLI
```
pip install awscli --upgrade --user
export PATH=~/.local/bin:$PATH (add into ~/.bashrc)
source ~/.bashrc
```
- Configure the CLI
```
aws configure
```
Then input the access key, secret access key, region time and output format.

- Download s3 objects
```
aws s3 sync s3://path filepath
```

