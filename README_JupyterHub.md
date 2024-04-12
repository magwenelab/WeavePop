## Reserve a VM
Go to `https://vcm.duke.edu/`
Click on "Reserve a VM"
Click on "Ubuntu Server 20.04"
Accept Terms and Conditions and wait to reveive an e-mail
Open a Terminal and access your VM by typing the SSH command in the e-mail. `ssh cz192@vcm-40396.vm.duke.edu`. If you have SSH Keys configured you don't need to provide your Duke password.
You will get this message:
```
The authenticity of host 'vcm-40528.vm.duke.edu (152.3.56.214)' can't be established.
ED25519 key fingerprint is SHA256:WKqxm4gl3ogU3vXzqCI5bUyujZkUwpO9o4FItfBSUXo.
This key is not known by any other names
Are you sure you want to continue connecting (yes/no/[fingerprint])?
```
Type `yes` and click Enter.
Now you are connected.
The number inside parentheses in the first line are the Public IP adress. We will need it to connect to JupyterHub.

## Install The Littelest JupyterHub

While being connected to the VM in your terminal in the admin user's home directory run:

```
sudo apt install python3 python3-dev git curl
```
You will be prompted for a password (use your Duke password) and then for a confirmation.

Write a text file called `requirements.txt` with this lines inside:
```
pandas
biopython
click
duckdb # SPECIFY SAME VERSION AS THE ONE USED TO BUILD THE DATABASE
```
Type the command bellow and replace the `<admin-user-name>` with the name you want your username to be and run it. It will take 5 to 10 minutes and it will say `Done!` when it finishes.  

```
curl -L https://tljh.jupyter.org/bootstrap.py | sudo -E python3 - --admin <admin-user-name> --user-requirements-txt requirements.txt 
```
From: `https://tljh.jupyter.org/en/latest/install/custom-server.html`

## Setup JupyterHub
In a browser go to `https://152.3.56.214`.  
Enter the admin username and create a password for it.
Using the button with the plus "+" symbol you will access the "Launcher". From there you can open a terminal, notebook, python console, etc.  

### Create users

In the graphic interface:
Go to "Files/Hub Control Panel".  
Go to "Admin", there you can add users.

### Upload data to the server

To upload files to the JupyterHub there are two options:
1. You can use the "Upload" button in the server to get files from your local machine to your Jupyter home. 
2. Using `sudo -E` (either from server's terminal or a local terminal connected with SSH to the VM) you can move files from the VM admin user's home `/home/admin-user-name` to wherever you need it (see below).


### Share data with all users

The way JupyterHub recommends to share data is by creating a folder in `/srv/data/my_shared_data_folder` and distribute it to users by making a symbolic link to it in `/etc/skel/`. Only users created after running this will see the folder in their home. Files moved to or created in that folder in any moment will be readable by all users that see the folder.  
```
sudo mkdir -p /srv/data/my_shared_data_folder
cd /etc/skel
sudo ln -s /srv/data/my_shared_data_folder my_shared_data_folder
```
The problem with that is that in our VM we don't have enough storage in `/srv`. So we can use `home/` to put the shared directory and do the symbolic link to it in `/etc/skel/`.  
To have the `Tutorial.md` and `Database_description.md` visible to all users directly in their home (without them having to look for them in the shared folder) they can be placed directly in `/home/` and link them individualy to `/etc/skel/`.  

Make sure the permissions are set to read-only for other users with:
```
sudo -E chmod 644 my_shared_data_folder/* Tutorial.md Database_description.md
```
Note: DuckDB connection to a database will need to specify `read_only = True`. Otherwise it will complain that it does not have permissions.

### Markdown Preview as default for all users

You can change the settings going to "Settings/Settings Editor". In the left you can see a panel with the things whose settings you can change and you can click on it and change it manually. Then if you click in the right corner to "JSON Settings Editor" you will see a JSON format version of the settings, including what you changed. From there you can copy the parts that you want to change in the settings of all users. 

Make the following directory and file:
```
sudo mkdir -p /opt/tljh/user/share/jupyter/lab/settings
sudo nano /opt/tljh/user/share/jupyter/lab/settings/overrides.json
```

To make the Markdown Preview the default viewer for markdown files, paste this snippet in the `overrides.json` file:

```
{
"@jupyterlab/docmanager-extension:plugin":
{
    "autosave": true,
    "autosaveInterval": 120,
    "confirmClosingDocument": false,
    "lastModifiedCheckMargin": 500,
    "defaultViewers": {
        "markdown": "Markdown Preview"
    },
    "renameUntitledFileOnSave": true
}
}
```
Reload the configuration.
```
sudo tljh-config reload
```
From: `https://tljh.jupyter.org/en/latest/howto/user-env/override-lab-settings.html`


### Setup HTTPS
`https://tljh.jupyter.org/en/latest/howto/admin/https.html`

## May be usefull in the future
`https://github.com/jupyterhub/the-littlest-jupyterhub/issues/429`
`https://blog.umd.edu/crstl/useful-links/how-to-setup-and-use-your-own-conda-environment-in-jupyterhub/`
