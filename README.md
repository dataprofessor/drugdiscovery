# drugdiscovery

[How to build a web app for Drug Discovery in Python | Streamlit #26](https://youtu.be/0rqIwSeUImo)

<a href="https://youtu.be/0rqIwSeUImo"><img src="http://img.youtube.com/vi/0rqIwSeUImo/0.jpg" alt="How to build a web app for Drug Discovery in Python | Streamlit #26" title="How to build a web app for Drug Discovery in Python | Streamlit #26" width="400" /></a>

# Reproducing this web app
To recreate this web app on your own computer, do the following.

### Create conda environment
Firstly, we will create a conda environment called *drugdiscovery*
```
conda create -n drugdiscovery python=3.7.9
```
Secondly, we will login to the *drugdiscovery* environement
```
conda activate drugdiscovery
```
### Install prerequisite libraries

Download requirements.txt file

```
wget https://raw.githubusercontent.com/dataprofessor/drugdiscovery/main/requirements.txt

```

Pip install libraries
```
pip install -r requirements.txt
```
###  Download and unzip contents from GitHub repo

Download and unzip contents from https://github.com/dataprofessor/drugdiscovery/archive/main.zip

###  Launch the app

```
streamlit run app.py
```
