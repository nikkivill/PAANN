document.addEventListener("DOMContentLoaded", function () {
    // Define the populations array
    const populations = ["BEB", "PJL"]; // Example population list

    // Tab switching functionality
    const tabItems = document.querySelectorAll('.tab-item');
    const tabContents = document.querySelectorAll('.tab-content');
    
    // variable to store chart instances
    const chartInstances = {};
    //  variable to track the highlighted point
    let highlightedPosition = null;
    
    tabItems.forEach(item => {
        item.addEventListener('click', function() {
            // Remove active class from all tabs and tab contents
            tabItems.forEach(tab => tab.classList.remove('active'));
            tabContents.forEach(content => content.classList.remove('active'));
            
            // Add active class to clicked tab and its content section
            this.classList.add('active');
            
            // Show corresponding content
            const tabId = this.getAttribute('data-tab');
            document.getElementById(tabId).classList.add('active');
            
            // If the visualization tab is activated, initialize the charts
            if (tabId === 'viz-tab') {
                initializeCharts();
            }
        });
    });

    // Add click event to table rows
    function addTableRowClickHandlers() {
        const tableRows = document.querySelectorAll('#snpTable tbody tr');

        tableRows.forEach(row => {
            row.addEventListener('click', function () {
                // Get position value from the clicked row (position is in the third column)
                const position = this.cells[2].innerText;

                // Remove highlighting from previously selected row
                document.querySelectorAll('#snpTable tbody tr.highlighted').forEach(r => {
                    r.classList.remove('highlighted');
                });

                // Highlight the clicked row
                this.classList.add('highlighted');

                // Check if the charts exist, and if not, initialize them
                if (Object.keys(chartInstances).length === 0) {
                    initializeCharts();  // Initialize charts before switching tabs

                    // Delay the highlighting slightly to ensure charts are ready
                    setTimeout(() => {
                        highlightChartPoints(position);
                    }, 500);
                } else {
                    highlightChartPoints(position);
                }

                // Activate the visualization tab
                document.querySelector('[data-tab="viz-tab"]').click();
            });
        });
    }

    // Function to highlight the selected position in all stats charts
    function highlightChartPoints(position) {
        highlightedPosition = position;
        
        // Updated all chart instances with highlighting
        Object.keys(chartInstances).forEach(chartId => {
            updateChartHighlighting(chartId, position);
        });
    }

    function updateChartHighlighting(chartId, position) {
        const chart = chartInstances[chartId];
        if (!chart) return;
        
        // Ensure the position is the same type as the labels
        const positionIndex = chart.data.labels.indexOf(position.toString()); // Convert position to string if necessary
        console.log(`Position: ${position}, Index: ${positionIndex}`); // Debugging log
        
        // Reset all point styles
        chart.data.datasets[0].pointBackgroundColor = chart.data.datasets[0].data.map(() => chart.data.datasets[0].borderColor);
        chart.data.datasets[0].pointBorderColor = chart.data.datasets[0].data.map(() => chart.data.datasets[0].borderColor);
        chart.data.datasets[0].pointRadius = chart.data.datasets[0].data.map(() => 2);
        chart.data.datasets[0].pointHoverRadius = chart.data.datasets[0].data.map(() => 4);
        
        // Highlight selected point
        if (positionIndex !== -1) {
            chart.data.datasets[0].pointBackgroundColor[positionIndex] = 'red'; // Highlight point in red
            chart.data.datasets[0].pointBorderColor[positionIndex] = 'red';
            chart.data.datasets[0].pointRadius[positionIndex] = 5; // Increase point size for visibility
            chart.data.datasets[0].pointHoverRadius[positionIndex] = 7;
        }
        
        chart.update();
    }

    // Initialize charts function
    function initializeCharts() {
        if (searchType === "genomic_region") {
            if (allResults) {
                console.log("Initializing charts with data:", allResults);
    
                function createChart(canvasId, label, data, title, color) {
                    console.log(`Creating chart ${canvasId} with data:`, data);
    
                    var chartInstance = Chart.getChart(canvasId);
                    if (chartInstance) {
                        chartInstance.destroy();
                    }
    
                    var ctx = document.getElementById(canvasId);
                    if (ctx) {
                        const newChart = new Chart(ctx, {
                            type: "line",
                            data: {
                                labels: snpPositions,
                                datasets: [
                                    {
                                        label: label,
                                        data: data,
                                        borderColor: color.replace("0.5", "1"),
                                        backgroundColor: color,
                                        borderWidth: 1,
                                        fill: false,
                                        pointRadius: 2,
                                        pointBackgroundColor: data.map(() => color.replace("0.5", "1")),
                                        pointBorderColor: data.map(() => color.replace("0.5", "1")),
                                        pointHoverRadius: 4
                                    }
                                ]
                            },
                            options: {
                                responsive: true,
                                maintainAspectRatio: false,
                                interaction: {
                                    mode: "index",
                                    intersect: false
                                },
                                plugins: {
                                    legend: {
                                        display: false
                                    },
                                    title: {
                                        display: false
                                    },
                                    tooltip: {
                                        enabled: true,
                                        callbacks: {
                                            label: function(tooltipItem) {
                                                var index = tooltipItem.dataIndex;
                                                var rsid = snpRSIDs[index];  
                                                var position = snpPositions[index];
                                                var value = data[index];
                                                return `RSID: ${rsid} | Position: ${position} | Value: ${value}`;
                                            }
                                        }
                                    },
                                    zoom: { 
                                        pan: { 
                                            enabled: true,
                                            mode: 'xy'
                                        },
                                        zoom: {
                                            wheel: {
                                                enabled: true,
                                            },
                                            pinch: {
                                                enabled: true
                                            }, 
                                            mode: 'xy',
                                            scaleMode: 'xy'
                                        } 
                                    }
                                },
                                scales: {
                                    x: {
                                        title: { 
                                            display: true,
                                            text: "SNP Position" 
                                        },
                                        ticks: { autoSkip: true, maxTicksLimit: 10 }
                                    },
                                    y: { 
                                        title: { display: true, text: title },
                                        beginAtZero: false
                                    }
                                }
                            }
                        });
                        
                        chartInstances[canvasId] = newChart;
                        
                        if (highlightedPosition !== null) {
                            updateChartHighlighting(canvasId, highlightedPosition);
                        }
                    }
                }
    
                if (populations.includes("BEB")) {
                    createChart("tajimaBEBChart", "Tajima's D (BEB)", tajimaBEB, "Value", "rgba(54, 162, 235, 0.5)");
                    createChart("nucleotideBEBChart", "Nucleotide Diversity (BEB)", nucleotideBEB, "Value", "rgba(153, 102, 255, 0.5)");
                    createChart("haplotypeBEBChart", "Haplotype Diversity (BEB)", haplotypeBEB, "Value", "rgba(201, 203, 207, 0.5)");
                }
    
                if (populations.includes("PJL")) {
                    createChart("tajimaPJLChart", "Tajima's D (PJL)", tajimaPJL, "Value", "rgba(255, 99, 132, 0.5)");
                    createChart("nucleotidePJLChart", "Nucleotide Diversity (PJL)", nucleotidePJL, "Value", "rgba(75, 192, 192, 0.5)");
                    createChart("haplotypePJLChart", "Haplotype Diversity (PJL)", haplotypePJL, "Value", "rgba(255, 159, 64, 0.5)");
                }
                
                createChart("fstChart", "FST Score", fstScores, "Value", "rgba(255, 206, 86, 0.5)");
                createChart("xpehhChart", "XPEHH Score", xpehhScores, "Value", "rgba(75, 192, 192, 0.5)");
                
                addTableRowClickHandlers();
            } else {
                console.log("No results available.");
            }
        }
    }

    // Function to download chart as a PNG image
    function downloadChart(chartId, fileName) {
        const chartCanvas = document.getElementById(chartId);
        if (chartCanvas) {
            const image = chartCanvas.toDataURL("image/png");
            const link = document.createElement("a");
            link.href = image;
            link.download = fileName;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
    }
    
    // Adds some styling for highlighted table rows
    const style = document.createElement('style');
    style.textContent = `
        #snpTable tbody tr.highlighted {
            background-color: #f0f8ff;
            cursor: pointer;
        }
        #snpTable tbody tr {
            cursor: pointer;
        }
        #snpTable tbody tr:hover {
            background-color: #f5f5f5;
        }
    `;
    document.head.appendChild(style); // Add the style to the document
    
    // Add download buttons to chart wrappers
    const chartWrappers = document.querySelectorAll('.chart-wrapper');
    chartWrappers.forEach(wrapper => {
        const canvas = wrapper.querySelector('canvas');
        if (canvas) {
            const chartId = canvas.id; // Get chart ID
            const chartTitle = wrapper.querySelector('.chart-title').textContent; // Get chart title
            const fileName = chartTitle.replace(/\s+/g, '_') + '.png'; // Replace spaces with underscores
            
            const downloadButton = document.createElement('button'); // Create download button
            downloadButton.className = 'btn btn-download'; // Add Bootstrap classes
            downloadButton.innerHTML = '<i class="bi bi-download"></i> Download'; // Add download icon
            downloadButton.onclick = function() {
                downloadChart(chartId, fileName); // Call download function
            };
            wrapper.appendChild(downloadButton); // Add button to chart wrapper
        }
    });
});