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
const pmnsElement = document.getElementById("pmns");

// Mixing Angle Experimental Values, degrees
const t12 = 33.44;
const t23 = 49;
const t13 = 8.57;

// squared mass differences, 1e-4 eV^2
const ms12 = .75;
const ms23 = 24.4;

// default settings
const defaultNeutrinoData =
      {
          t12,
          t23,
          t13,
          ms12,
          ms23,
          L_max: 100,
          is_e: 1,
          is_m: 0,
          is_t: 0,
          E: 6 // MeV
      };


const control_categories = {
    angles: {
        title: 'Mixing Angles',
    },

    masses: {
        title: 'Squared Mass Differences',
    },

    state: {
        title: 'Initial State',
    },

    observatory: {
        title: 'Observatory Settings',
    },
};

const neutrino_controls = {
    t12: {
        title: '\\(\\theta_{12}\\)',
        unit: '\\(^\\circ\\)',
        category: 'angles',
        control: {
            type: 'slider',
            min: 0,
            max: 90,
            step: .1,
        }
    },

    t13: {
        title: '\\(\\theta_{13}\\)',
        unit: '\\(^\\circ\\)',
        category: 'angles',
        control: {
            type: 'slider',
            min: 0,
            max: 90,
            step: .1,
        }
    },

    t23: {
        title: '\\(\\theta_{23}\\)',
        unit: '\\(^\\circ\\)',
        category: 'angles',
        control: {
            type: 'slider',
            min: 0,
            max: 90,
            step: .1,
        }
    },

    ms12: {
        title: '\\(\\Delta m_{12}^2\\)',
        unit: '\\(\\cdot 10^{-4}\\text{eV}^2\\)',
        category: 'masses',
        control: {
            type: 'slider',
            min: 0,
            max: 50,
            step: .01,
        }
    },

    ms23: {
        title: '\\(\\Delta m_{23}^2\\)',
        unit: '\\(\\cdot 10^{-4}\\text{eV}^2\\)',
        category: 'masses',
        control: {
            type: 'slider',
            min: 0,
            max: 50,
            step: .01,
        }
    },

    L_max: {
        title: '\\(L_\\text{max}\\)',
        unit: ' \\(\\text{km}\\)',
        category: 'observatory',
        control: {
            type: 'slider',
            min: 0,
            max: 1000,
            step: .1,
        }
    },

    E: {
        title: '\\(E\\)',
        unit: ' \\(\\text{MeV}\\)',
        category: 'observatory',
        control: {
            type: 'slider',
            min: 0,
            max: 100,
            step: .1,
        }
    },

    is_e: {
        title: '\\(e\\)',
        unit: '',
        category: 'state',
        control: {
            type: 'slider',
            min: 0,
            max: 1,
            step: .01,
        }
    },

    is_m: {
        title: '\\(\\mu\\)',
        unit: '',
        category: 'state',
        control: {
            type: 'slider',
            min: 0,
            max: 1,
            step: .01,
        }
    },

    is_t: {
        title: '\\(\\tau\\)',
        unit: '',
        category: 'state',
        control: {
            type: 'slider',
            min: 0,
            max: 1,
            step: .01,
        }
    },
}

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

function propagation_step(pmns_matrix, ms12, ms23, dL, E) {
    E = E / 1000;
    ms12 = ms12 * 1e-4;
    ms23 = ms23 * 1e-4;

    const propagator = math.matrix([[math.exp(math.complex(0, -ms12 * 2.54 * dL / E)), 0, 0],
                                    [0, 1, 0],
                                    [0, 0, math.exp(math.complex(0, ms23 * 2.54 * dL / E))]]);

    return math.multiply(pmns_matrix, propagator, math.transpose(pmns_matrix));
}

function propagate(initial_state, pmns_matrix, ms12, ms23, L, E) {
    const propagator = propagation_step(pmns_matrix, ms12, ms23, L, E);

    return math.multiply(propagation_step, initial_state);
}

function state_to_probabilies(state) {
    return state.map(value => math.multiply(value.conjugate(), value).re);
}

function plot_propagation(neutrino_data) {
    const {t12, t23, t13, ms12, ms23, L_max, E, is_e, is_m, is_t} = neutrino_data;

    const norm = math.sqrt(is_e * is_e + is_m * is_m + is_t * is_t);
    let state = math.matrix([
        [is_e / norm],
        [is_m / norm],
        [is_t / norm]
    ]);

    const pmns_matrix = pmns(t12, t23, t13);

    const dL = ((L_max) / 1000);
    const lengths = math.range(0, L_max, dL);

    const common_options = {
        pointRadius: 0,
    };

    const datasets = [
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

    const propagator = propagation_step(pmns_matrix, ms12, ms23, dL, E)

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

///////////////////////////////////////////////////////////////////////////////
//                                  UI Stuff                                 //
///////////////////////////////////////////////////////////////////////////////

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


    pmnsElement.textContent = result;
    renderMathInElement(pmnsElement);
}

function setUpApp() {
    const neutrino_data = {...defaultNeutrinoData, ...getWindowNeutrinoData()};
    renderControls(neutrino_data);

    const {t12, t23, t13} = getNeutrinoData();
    renderPmnsMatrix(pmns(t12, t23, t13));

    plot_propagation(neutrino_data);
}

document.addEventListener("DOMContentLoaded", setUpApp);
window.addEventListener("popstate", setUpApp);

function makeControlElement(id, element, default_value) {
    const {min, max, step} = element.control;

    const input_div = document.createElement("div");
    const label = document.createElement("label");
    const input = document.createElement("input");


    label.innerHTML = `${element.title} = <span id="val_${id}">${default_value}</span>${element.unit}`;
    label.htmlFor = id

    input.type = "range";
    input.step = step;
    input.min = min;
    input.max = max;
    input.id = id;
    input.value = default_value;
    input.label = label;

    input_div.className = "input-group vertical";
    input_div.appendChild(label);
    input_div.appendChild(input);

    return input_div;
}

function renderControls(neutrino_data) {
    controls.innerHTML = "";
    for (let category in control_categories) {
        const category_element = document.createElement("fieldset");
        category_element.id = category;

        const legend = document.createElement("legend");
        legend.textContent = control_categories[category].title;

        category_element.appendChild(legend);
        controls.appendChild(category_element);
    }

    for (let id in neutrino_controls) {
        const element = neutrino_controls[id];
        const group = controls.querySelector("#" + element.category);
        const control = makeControlElement(id, element, neutrino_data[id]);

        group.appendChild(control);
    }

    renderMathInElement(controls);
    controls.addEventListener('input', handleInput);
}

function getNeutrinoData() {
    const data = {};

    for (let element in neutrino_controls) {
        data[element] = parseFloat(controls.elements[element].value);
    }

    return {...defaultNeutrinoData, ...data};
}

function updateControlLabel(control) {
    control.label.childNodes.item(3).textContent = control.value;
}


function serializeNeutrinoData(data) {
    return encodeURI(JSON.stringify(data));
}

function getWindowNeutrinoData() {
    const data_string = decodeURI(window.location.href);
    return JSON.parse(data_string.slice(data_string.indexOf("?settings=") + "?settings=".length));
}

let historyTimeout = null;

function handleInput(event) {
    const control = event.target;
    normalizeSliders()

    const neutrino_data = getNeutrinoData();
    const current_href = window.location.href;

    window.clearTimeout(historyTimeout);
    historyTimeout = window.setTimeout(() => {
        window.clearTimeout(historyTimeout);
        window.history.pushState(neutrino_data, "",
                                 current_href.slice(0, current_href.indexOf("?"))
                                 + "?settings=" + serializeNeutrinoData(neutrino_data))
    }, 1000);

    const {t12, t23, t13} = neutrino_data;
    renderPmnsMatrix(pmns(t12, t23, t13));

    updateControlLabel(control);
    plot_propagation(neutrino_data);
}

function normalizeSliders() {
    const {is_e, is_m, is_t} = getNeutrinoData();
    const norm = math.sqrt(is_e * is_e + is_m * is_m + is_t * is_t);

    for (let name of ['is_e', 'is_t', 'is_m']) {
        const control = controls[name];
        control.value = control.value / norm
        updateControlLabel(control);
    }
}
