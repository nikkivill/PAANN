document.addEventListener("DOMContentLoaded", function () {
    // Makes sure run only after the page is loaded
    const populations = ["BEB", "PJL"]; // Example population list

    // Tab switching (selects all with .tab-item and .tab-content classes)
    const tabItems = document.querySelectorAll('.tab-item');
    const tabContents = document.querySelectorAll('.tab-content');

    // Stores chart instances
    const chartInstances = {};
    let highlightedPosition = null; // Tracks the highlighted position in the table
    let isChartReady = false; // Tracks if charts are initialised

    // Adds click event listener to each tab
    tabItems.forEach(item => {
        item.addEventListener('click', function () {
            // Removes the active class from all tabs and tab contents as a reset
            tabItems.forEach(tab => tab.classList.remove('active'));
            tabContents.forEach(content => content.classList.remove('active'));

            // Activates the clicked tab
            this.classList.add('active');

            // Activates the tab content with the same id as the clicked tab
            const tabId = this.getAttribute('data-tab');
            document.getElementById(tabId).classList.add('active');

            // If visualisation tab is activated, initialises the charts
            if (tabId === 'vis-tab') {
                initializeCharts();
            }
        });
    });

    // Adds click event listener to the table rows
    function addTableRowClickHandlers() {
        const tableRows = document.querySelectorAll('#snpTable tbody tr'); // Selects table rows, excluding header
        const tablePositions = Array.from(tableRows).map(row => row.cells[2].innerText.trim()); // Extracts positions from each row
        console.log("Table Positions:", tablePositions);

        tableRows.forEach(row => {
            row.addEventListener('click', function () {
                const position = this.cells[2].innerText.trim(); // Gets the position of the clicked row
                console.log("Clicked Position:", position);

                // Removes .highlighted class from previously selected rows
                document.querySelectorAll('#snpTable tbody tr.highlighted').forEach(r => r.classList.remove('highlighted'));

                // Adds .highlighted class to the clicked row
                this.classList.add('highlighted');

                // Initialises charts if not ready
                if (!isChartReady) {
                    console.log("Charts not ready yet. Initializing...");
                    initializeCharts();
                    setTimeout(() => highlightChartPoints(position), 500); // Delay highlighting to ensure charts are ready
                } else {
                    highlightChartPoints(position);
                }

                // Activates the visualisation tab
                document.querySelector('[data-tab="vis-tab"]').click();
            });
        });
    }

    // Highlights the selected position in all stats charts
    function highlightChartPoints(position) {
        highlightedPosition = position;
        Object.keys(chartInstances).forEach(chartId => updateChartHighlighting(chartId, position));
    }

    // Updates chart highlighting for a specific position
    function updateChartHighlighting(chartId, position) {
        const chart = chartInstances[chartId];
        if (!chart) {
            console.error(`Chart ${chartId} not found!`);
            return;
        }
        
        // Turns SNP position into string, trims and indexes (was issue before as if not done, string integer mismatch occurs)
        const normalizedPosition = position.toString().trim();
        const positionIndex = chart.data.labels.map(String).indexOf(normalizedPosition);
        console.log(`Position: ${normalizedPosition}, Index: ${positionIndex}, Labels:`, chart.data.labels);

        // Resets all point styles to default
        chart.data.datasets[0].pointBackgroundColor = chart.data.datasets[0].data.map(() => chart.data.datasets[0].borderColor);
        chart.data.datasets[0].pointBorderColor = chart.data.datasets[0].data.map(() => chart.data.datasets[0].borderColor);
        chart.data.datasets[0].pointRadius = chart.data.datasets[0].data.map(() => 2);
        chart.data.datasets[0].pointHoverRadius = chart.data.datasets[0].data.map(() => 4);

        // Highlights chosen SNP in red on charts
        chart.data.datasets[0].pointBackgroundColor[positionIndex] = 'red';
        chart.data.datasets[0].pointBorderColor[positionIndex] = 'red';
        chart.data.datasets[0].pointRadius[positionIndex] = 5;
        chart.data.datasets[0].pointHoverRadius[positionIndex] = 7;

        chart.update();
    }

    // Initialises charts if search type is genomic region
    function initializeCharts() {
        if (searchType === "genomic_region" && allResults) {
            console.log("Initializing charts with data:", allResults);
            console.log("Chart Labels (snpPositions):", snpPositions);

            // Helper function to create a chart
            function createChart(canvasId, label, data, title, color) {
                console.log(`Creating chart ${canvasId} with data:`, data);

                var chartInstance = Chart.getChart(canvasId);
                if (chartInstance) chartInstance.destroy();

                var ctx = document.getElementById(canvasId);
                if (ctx) {
                    const newChart = new Chart(ctx, {
                        type: "line",
                        data: {
                            labels: snpPositions,
                            datasets: [{
                                label: label,
                                data: data,
                                // Styling
                                borderColor: color.replace("0.5", "1"),
                                backgroundColor: color,
                                borderWidth: 1,
                                fill: false,
                                pointRadius: 2,
                                pointBackgroundColor: data.map(() => color.replace("0.5", "1")),
                                pointBorderColor: data.map(() => color.replace("0.5", "1")),
                                pointHoverRadius: 4
                            }]
                        },
                        // Makes the chart responsive and interactive
                        options: {
                            responsive: true,
                            maintainAspectRatio: false,
                            interaction: { mode: "index", intersect: false },
                            plugins: {
                                legend: { display: false },
                                title: { display: false },
                                tooltip: {
                                    enabled: true,
                                    callbacks: { // Custom tooltip label when hovering
                                        label: function (tooltipItem) {
                                            var index = tooltipItem.dataIndex;
                                            var rsid = snpRSIDs[index];
                                            var position = snpPositions[index];
                                            var value = data[index];
                                            return `RSID: ${rsid} | Position: ${position} | Value: ${value}`;
                                        }
                                    }
                                },
                                // Zoom and pan plugin for better chart navigation
                                zoom: {
                                    pan: { enabled: true, mode: 'xy' },
                                    zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: 'xy', scaleMode: 'xy' }
                                }
                            }, // Custom axes titles and ticks
                            scales: {
                                x: { title: { display: true, text: "SNP Position" }, ticks: { autoSkip: true, maxTicksLimit: 10 } },
                                y: { title: { display: true, text: title }, beginAtZero: false }
                            }
                        }
                    });

                    // Stores chart instance for later use
                    chartInstances[canvasId] = newChart;
                    if (highlightedPosition !== null) updateChartHighlighting(canvasId, highlightedPosition);
                }
            }

            // Creates charts for populations
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

            // Creates comparative stats charts
            createChart("fstChart", "FST Score", fstScores, "Value", "rgba(255, 206, 86, 0.5)");
            createChart("xpehhChart", "XPEHH Score", xpehhScores, "Value", "rgba(75, 192, 192, 0.5)");

            isChartReady = true; // Mark charts as ready
            addTableRowClickHandlers();
        } else {
            console.log("No results available.");
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

    // Adds styling for highlighted table rows, css
    const style = document.createElement('style');
    style.textContent = `
        #snpTable tbody tr.highlighted { background-color: #f0f8ff; cursor: pointer; }
        #snpTable tbody tr { cursor: pointer; }
        #snpTable tbody tr:hover { background-color: #f5f5f5; }
    `;
    document.head.appendChild(style);

    // Add download buttons to chart wrappers
    const chartWrappers = document.querySelectorAll('.chart-wrapper');
    chartWrappers.forEach(wrapper => { // iterates through each chart wrapper
        const canvas = wrapper.querySelector('canvas'); 
        if (canvas) {  // finds canvas element, if inside, gets ID, title (for named file) and adds <button> element with bootstrap styling
            const chartId = canvas.id;
            const chartTitle = wrapper.querySelector('.chart-title').textContent;
            const fileName = chartTitle.replace(/\s+/g, '_') + '.png';

            const downloadButton = document.createElement('button');  // Creates download button
            downloadButton.className = 'btn btn-download'; // Bootstrap styling
            downloadButton.innerHTML = '<i class="bi bi-download"></i> Download'; // Adds download icon
            downloadButton.onclick = () => downloadChart(chartId, fileName);
            wrapper.appendChild(downloadButton);// Add button to the chart wrapper
        }
    });
});