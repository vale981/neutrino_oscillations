///////////////////////////////////////////////////////////////////////////////
//                               Setup/Globals                               //
///////////////////////////////////////////////////////////////////////////////

const ctx = document.getElementById('neutrinorender').getContext('2d');
const render = new Chart(ctx, {
    type: 'line',
    options: {
        responsive: true,
        maintainAspectRatio: false,
        tooltips: {
            mode: 'index',
            intersect: false,
        },
        animation: {
            duration: 0,
        },
        scales: {
            x: {
                scaleLabel: {
                    display: true,
                    labelString: "Distance [km]"
                }
            },
            y: {
                scaleLabel: {
                    display: true,
                    labelString: "Survival Probability"
                }
            }
        }
    },
});

const controls = document.forms.controls;
const settings = document.getElementById('settings');

// Mixing Angle Experimental Values, degrees
const t12 = 33.44;
const t23 = 49;
const t13 = 8.57;

// squared mass differences, 1e-4 eV^2
const ms13 = 25.1;
const ms23 = 24.4;

// default settings
const defaultNeutrinoData =
      {
          t12,
          t23,
          t13,
          ms13,
          ms23,
          L_range: [0, 100],
          E: 6 // MeV
      };

///////////////////////////////////////////////////////////////////////////////
//                                    Util                                   //
///////////////////////////////////////////////////////////////////////////////
function deg_to_rad(x) {
    return (x / 180) * math.pi
}

///////////////////////////////////////////////////////////////////////////////
//                         Neutrino Mixing Specifics                         //
///////////////////////////////////////////////////////////////////////////////

/**
 * Compute the PMNS Matrix.
 *
 * Angle Arguments in Degrees.
 */
function pmns(t12, t23, t13) {
    t12 = deg_to_rad(t12);
    t13 = deg_to_rad(t13);
    t23 = deg_to_rad(t23);

    // Compute all sines, cosines
    let c12 = math.cos(t12);
    let c13 = math.cos(t13);
    let c23 = math.cos(t23);

    let s12 = math.sin(t12);
    let s13 = math.sin(t13);
    let s23 = math.sin(t23);

    return math.matrix([
        [c12 * c13, s12 * c13, s13],
        [-s12 * c23 - c12 * s23 * s13, c12 * c23 - s12 * s23 * s13, s23 * c13],
        [s12 * s23 - c12 * c23 * s13, -c12 * s23 - s12 * c23 * s13, c23 * c13]
    ]);
}

function propagation_step(pmns_matrix, ms13, ms23, dL, E) {
    E = E / 1000;
    ms13 = ms13 * 1e-4;
    ms23 = ms23 * 1e-4;

    const propagator = math.matrix([[math.exp(math.complex(0, -ms13 * 2.54 * dL / E)), 0, 0],
                                    [0, math.exp(math.complex(0, -ms23 * 2.54 * dL / E)), 0],
                                    [0, 0, 1]]);

    return math.multiply(pmns_matrix, propagator, math.transpose(pmns_matrix));
}

function propagate(initial_state, pmns_matrix, ms13, ms23, L, E) {
    const propagator = propagation_step(pmns_matrix, ms13, ms23, L, E);

    return math.multiply(propagation_step, initial_state);
}

function state_to_probabilies(state) {
    return state.map(value => math.multiply(value.conjugate(), value).re);
}

function plot_propagation(state, neutrino_data) {
    const {t12, t23, t13, ms13, ms23, L_range, E} = neutrino_data;

    pmns_matrix = pmns(t12, t23, t13);

    dL = ((L_range[1] - L_range[0]) / 1000);
    lengths = math.range(L_range[0], L_range[1], dL);

    common_options = {
        pointRadius: 0,
    };

    datasets = [
        {
            data: [],
            borderColor: 'green',
            label: 'electron',
            ...common_options
        },
        {
            data: [],
            borderColor: 'blue',
            label: 'muon',
            ...common_options
        },
        {
            data: [],
            borderColor: 'red',
            label: 'tauon',
            ...common_options
        }
    ];

    propagator = propagation_step(pmns_matrix, ms13, ms23, dL, E)

    for (const length of lengths._data) {
        state = math.multiply(propagator, state);

        probs = state_to_probabilies(state);
        probs.forEach((p, index) => {
            datasets[index[0]].data.push(p);
        })
    }

    render.data.labels = lengths.map(l => math.round(l, 2))._data;
    render.data.datasets = datasets;
    render.update({
        duration: 0,
        easing: 'easeInOutBack'
    });
}

function initForm() {
    controls.elements.t13.value = t13;
    controls.elements.t12.value = t12;
    controls.elements.t23.value = t23;

    controls.elements.L_max.value = defaultNeutrinoData.L_range[1];
    controls.elements.energy.value = defaultNeutrinoData.E;

    controls.elements.ms13.value = ms13;
    controls.elements.ms23.value = ms23;
}

function getNeutrinoData() {
    return {
        ...defaultNeutrinoData,
        t13: controls.elements.t13.value,
        t12: controls.elements.t12.value,
        t23: controls.elements.t23.value,
        ms13: controls.elements.ms13.value,
        ms23: controls.elements.ms23.value,
        E: controls.elements.energy.value,
        L_range: [0, controls.elements.L_max.value],
    };
}

function renderPmnsMatrix(pmns_matrix) {
    let result = "$$U=\\begin{pmatrix}\n";

    for (let row of pmns_matrix._data) {
        for (let val of row) {
            result += `${math.round(val, 3)} &`
        }
        result = result.slice(0, -1);
        result += "\\\\\n";
    }

    result += "\\end{pmatrix}$$";

    return result;
}

function displayCurrentSettings() {
    const {t12, t23, t13, ms13, ms23, L_range, E} = getNeutrinoData();
    pmns_matrix = pmns(t12, t23, t13);

    // render pmns
    settings.querySelector("#pmns").textContent = renderPmnsMatrix(pmns_matrix);

    settings.querySelector("#sliders_1").innerHTML = `
\\(\\theta_{12} = ${t12}^\\circ\\)
<br>
\\(\\theta_{13} = ${t13}^\\circ\\)
<br>
\\(\\theta_{23} = ${t23}^\\circ\\)
    `;
    settings.querySelector("#sliders_2").innerHTML = `
\\(\\Delta m_{13}^2 = ${ms13}\\cdot 10^{-4}\\text{GeV}\\)
<br>
\\(\\Delta m_{23}^2 = ${ms23}\\cdot 10^{-4}\\text{GeV}\\)
<br>
\\(L_\\text{max} = ${L_range[1]}\\text{km}\\)
<br>
\\(E = ${E}\\text{MeV}\\)
    `;

    renderMathInElement(settings);

};

document.addEventListener("DOMContentLoaded", (event) => {
    initForm();
    displayCurrentSettings();
    plot_propagation(math.matrix([
        [1],
        [0],
        [0]
    ]), defaultNeutrinoData);
});

controls.addEventListener('input', () => {
    displayCurrentSettings();
    plot_propagation(math.matrix([
        [1],
        [0],
        [0]
    ]), getNeutrinoData());
});
