let buffer = {};
let config = {};
let next = 0;

document.addEventListener("DOMContentLoaded", () => {
    const socket = io();
    let progressChart;
    const progressData = {
        labels: [],
        datasets: [
            {
                label: 'CAI',
                borderColor: 'rgba(75, 192, 192, 1)',
                backgroundColor: 'rgba(75, 192, 192, 0.2)',
                data: [],
                yAxisID: 'y-left' // Attach CAI to the left Y-axis
                // hidden: true // Hide initially
            },
            {
                label: 'AUP',
                borderColor: 'rgba(255, 99, 132, 1)',
                backgroundColor: 'rgba(255, 99, 132, 0.2)',
                data: [],
                yAxisID: 'y-right', // Attach AUP to the right Y-axis
                hidden: true // Hide initially
            },
            {
                label: 'EFE',
                borderColor: 'rgba(54, 162, 235, 1)',
                backgroundColor: 'rgba(54, 162, 235, 0.2)',
                data: [],
                yAxisID: 'y-right', // Attach AUP to the right Y-axis
                hidden: true // Hide initially
            }
        ]
    };

    const ctx = document.getElementById('progressChart').getContext('2d');
    progressChart = new Chart(ctx, {
        type: 'line',
        data: progressData,
        options: {
            maintainAspectRatio: false,
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Step'
                    }
                },
                'y-left': {
                    type: 'linear',
                    position: 'left',
                    title: {
                        display: true,
                        text: 'CAI'
                    },
                    ticks: {
                        beginAtZero: true,
                        max: 1 // Adjust this based on your CAI data range
                    }
                },
                'y-right': {
                    type: 'linear',
                    position: 'right',
                    title: {
                        display: true,
                        text: 'AUP / EFE'
                    },
                    ticks: {
                        beginAtZero: true,
                        max: 1 // Adjust this based on your AUP/EFE data range
                    },
                    grid: {
                        drawOnChartArea: false // Only gridlines for the left Y-axis
                    }
                }
            },
            plugins: {
                annotation: {
                    annotations: {
                        line1: {
                            type: 'line',
                            yMin: 0.8,
                            yMax: 0.8,
                            borderColor: 'rgba(255, 99, 132, 0.5)',
                            borderWidth: 2,
                            scaleID: 'y-left',
                            label: {
                                content: 'CAI Threshold (0.8)',
                                enabled: true,
                                position: 'end'
                            }
                        }
                        // Add more horizontal lines here if needed
                    }
                }
            }
        }
    });

    let next = 0;

    socket.on('connect', () => {
        console.log('Connected to server');
    });

    socket.on('task_update', (response) => {
        // console.log(response);
        buffer[response.step] = response;
        next = processBuffer(buffer, next);
    });

    function processBuffer(buffer, next) {
        while (next in buffer) {
            applyUpdate(buffer[next]);
            delete buffer[next];
            next++;
        }
        return next;
    }

    function applyUpdate(response) {
        updateChart(response);
        updateCDS(response.cds);
        updateStep(response.step);
    }

    socket.on('arwa_sync_error', (response) => {
        const data = response.result;
        console.error(data.error);
    });

    function updateChart(response) { 
        const {
            step,
            measures,
        } = response;
        progressData.labels.push(step);
        ["CAI", "AUP", "EFE"].forEach((key, i) => {
            if (key in measures) {
                progressData.datasets[i].data.push(measures[key]);
            }
        });
        progressChart.update();
    }
    function updateCDS(cds) {
        const cdsContainer = document.getElementById('cds');
        cdsContainer.innerHTML = cds.join(', ');
    }

    function updateStep(step) {
        const currentStepElement = document.getElementById('currentStep');
        currentStepElement.textContent = step;
    }

    function resetChart(config) {
        let {
            steps,
            stability,
            cai_threshold,
            cai_exp_scale,
        } = config
        // Show relevant datasets based on selection
        progressData.datasets[0].hidden = false; // CAI is always shown
        if (stability === 'aup') {
            progressData.datasets[1].hidden = false;
            progressData.datasets[2].hidden = true;
            progressChart.options.scales['y-right'].title.text = 'AUP';
            progressChart.options.scales['y-right'].display = true;
        } else if (stability === 'efe') {
            progressData.datasets[1].hidden = true;
            progressData.datasets[2].hidden = false;
            progressChart.options.scales['y-right'].title.text = 'EFE';
            progressChart.options.scales['y-right'].display = true;
        } else {
            progressData.datasets[1].hidden = true;
            progressData.datasets[2].hidden = true;
            // remove the right Y-axis, or hide it
            progressChart.options.scales['y-right'].display = false;
        }
        progressData.labels.length = 0;
        progressData.datasets.forEach(dataset => {
            dataset.data.length = 0;
            // dataset.hidden = true; // Hide all initially
        });
        progressChart.update();
    }
    
    document.getElementById('arwaForm').addEventListener('submit', (event) => {
        event.preventDefault();
        updateStep(0);
        next = 0;
        const formData = new FormData(event.target);
        const requestData = Object.fromEntries(formData.entries());
        requestData.cai_threshold = parseFloat(requestData.cai_threshold);
        requestData.cai_exp_scale = parseFloat(requestData.cai_exp_scale);
        requestData.steps = parseInt(requestData.steps, 10);
        config = {
            steps: requestData.steps,
            stability: requestData.stability,
            cai_threshold: requestData.cai_threshold,
            cai_exp_scale: requestData.cai_exp_scale,
        }
        resetChart(config);
        socket.emit('arwa_websocket_celery', requestData);
    });

    document.getElementById('test').addEventListener('click', () => {
        document.getElementById('cai_threshold').value = 0.9;
        document.getElementById('cai_exp_scale').value = 1.0;
        document.getElementById('freq_table_path').value = "homosapiens.txt";
        document.getElementById('stability').value = "efe";
        document.getElementById('aa_seq').value = "MVSKGEELFTGVVPILVELDGDVNGH";
        document.getElementById('steps').value = 10;
    });
});

function openPopup() {
    document.getElementById('popup').style.display = 'block';
    document.getElementById('overlay').style.display = 'block';
}

function closePopup() {
    document.getElementById('popup').style.display = 'none';
    document.getElementById('overlay').style.display = 'none';
}

function searchPopup() {
    // Define the taxID you want to use
    // const taxID = '12345'; // Replace with the actual taxID you want to pass
    const taxID = document.getElementById('species').value;
    // Construct the URL with the taxID
    const url = `/custom_codon_table/${taxID}`;

    // Use fetch to make a GET request to the route
    fetch(url, {
        method: 'GET',
    })
    .then(response => {
        if (!response.ok) {
            alert('No data found 1');
            throw new Error('Network response was not ok ' + response.statusText);
        }
        return response.text(); // Or response.json() if the response is JSON
    })
    .then(data => {
        console.log(data); // Handle the response data
        // You can update your UI or do something else with the data here
        if (data) {
            // Check if taxID exists in the options
            let optionExists = false;
            var select = document.getElementById('freq_table_path');
            for (var i = 0; i < select.options.length; i++) {
                if (select.options[i].value === taxID) {
                    optionExists = true;
                    break;
                }
            }

            // If the option doesn't exist, create a new option
            if (!optionExists) {
                var newOption = new Option(taxID, taxID);
                select.add(newOption);
            }

            // Set the value of the select element to taxID
            select.value = taxID;
            document.getElementById('freq_table_path').value = taxID;
            closePopup();
        } else {
            alert('No data found');
        }
    })
    .catch(error => {
        console.error('There was a problem with the fetch operation:', error);
    });
}