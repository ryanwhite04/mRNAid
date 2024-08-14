---
numbersections: true
figPrefix: "Figure "
fig-suffix: ": "
...

# Admin Guide for mRNAid Web App Deployment

This manual will guide you through the setup and deployment of the mRNAid web application, covering both development and production environments.

## Prerequisites

1. **Install Node.js and npm:**

    Nodejs needs to be installed to run the react server

    Installation instructions are [here](https://nodejs.org/en/learn/getting-started/how-to-install-nodejs)

1. **Install Conda or Miniconda**

    Conda is required to run the flask server, and task queue

    Installation instructions for your specific system can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

1. **Clone the Repository**

    The repository is publically available and can be found [here](https://github.com/ryanwhite04/mRNAid)

1. **Create the conda environment**

    ```sh
    conda create --name mRNAid python3.9    
    ```
1. **Install python requirements**

    Install all the python packages in requiremnts.txt using conda

1. **Install redis**

    To install redis, cd into /mRNAid/mrnaid/backend/flask_app and run

    ```sh
    make redis-install
    ```

1. **Install npm packages**

    navigate to /mRNAid/mrnaid/frontend/

    ```sh
    npm install
    ```

1. **Start the servers**

    navigate back to /mRNAid/mrnaid/backend/flask_app

    and start each server with the corresponding make command

    ```sh
    make redis-run # to start redis server
    make celery-run # to spawn celery workers
    make flower-run # to launch flower UI on port 5566
    make flask-run # to start backend server on port 5000
    make react-run # to launch react server on port 3000
    ```

    Each of the servers must be started in it's own terminal
    They must also be started in the specific order outlined above

1. **Server Diagram**

    A diagram of the servers and there interaction is shown in Figure 1 below

    ![Servers](images/servers.png)

1. **Submit an optimization task**

    Navigate the [local site]("http://localhost:3000) to view the frontend shown in figure 2

    ![React Frontend Empty](images/react_frontend_empty.png)

    You can hit example to fill in an example request for testing shown in figure 3

    ![React frontend full](images/react_frontend_full.png)

    Hitting submit will take you to the loading page shown in figure 4

    ![React frontend loading](images/react_frontend_loading.png)

    And after a few second the results will appear as in figure 5

    ![React frontend finished](images/react_frontend_finished.png)

    You can also open the [flower frontend]("http://localhost:5566) shown in figure 6

    ![Flower frontend](images/flower_frontend.png)  

1. **Appendix**