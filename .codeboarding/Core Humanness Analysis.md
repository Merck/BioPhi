```mermaid
graph LR
    User_Interface_CLI["User Interface & CLI"]
    Task_Orchestration["Task Orchestration"]
    Humanness_Data_Models["Humanness Data Models"]
    Core_Humanness_Analysis["Core Humanness Analysis"]
    Humanization_Core["Humanization Core"]
    User_Interface_CLI -- "initiates" --> Task_Orchestration
    Task_Orchestration -- "orchestrates" --> Core_Humanness_Analysis
    Task_Orchestration -- "orchestrates" --> Humanization_Core
    Task_Orchestration -- "processes" --> Humanness_Data_Models
    Core_Humanness_Analysis -- "calculates using" --> Humanness_Data_Models
    Humanization_Core -- "processes/produces" --> Humanness_Data_Models
    Humanness_Data_Models -- "structures data for" --> Core_Humanness_Analysis
    Humanness_Data_Models -- "structures data for" --> Humanization_Core
```
[![CodeBoarding](https://img.shields.io/badge/Generated%20by-CodeBoarding-9cf?style=flat-square)](https://github.com/CodeBoarding/GeneratedOnBoardings)[![Demo](https://img.shields.io/badge/Try%20our-Demo-blue?style=flat-square)](https://www.codeboarding.org/demo)[![Contact](https://img.shields.io/badge/Contact%20us%20-%20contact@codeboarding.org-lightgrey?style=flat-square)](mailto:contact@codeboarding.org)

## Component Details

This architecture outlines the BioPhi humanization subsystem, which processes antibody sequences to assess and improve their 'humanness'. The `User Interface & CLI` component serves as the entry point, receiving user input. This input is then passed to the `Task Orchestration` component, which manages the asynchronous execution of either humanness analysis or humanization tasks. The `Core Humanness Analysis` component performs the detailed calculation of humanness scores using methods like OASis, relying on data structures defined in `Humanness Data Models`. Similarly, the `Humanization Core` component modifies antibody sequences to increase humanness, also interacting with `Humanness Data Models` for data representation. The `Humanness Data Models` component acts as a central data repository, structuring information for both analysis and humanization processes.

### User Interface & CLI
This component handles user interactions, including receiving antibody sequences and humanization parameters via web forms or command-line interfaces. It is responsible for initiating the core humanization and humanness calculation processes.


**Related Classes/Methods**:

- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/cli/oasis.py#L19-L96" target="_blank" rel="noopener noreferrer">`biophi.humanization.cli.oasis:oasis` (19:96)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/cli/sapiens.py#L33-L124" target="_blank" rel="noopener noreferrer">`biophi.humanization.cli.sapiens:sapiens` (33:124)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/views.py#L97-L160" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.views:humanize_post` (97:160)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/views.py#L333-L361" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.views:humanness_post` (333:361)</a>


### Task Orchestration
This component manages the asynchronous execution of humanization and humanness calculation tasks. It schedules tasks, monitors their progress, and provides methods for retrieving and formatting results for display or export.


**Related Classes/Methods**:

- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L32-L43" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks.HumanizeAntibodyTaskResult:to_series` (32:43)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L70-L89" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks.HumanizeAntibodyTaskResult:to_sheets` (70:89)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L123-L171" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks:humanize_antibody_task` (123:171)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L175-L205" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks:mutate_humanized_antibody_task` (175:205)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L220-L244" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks.HumannessTaskResult:to_series` (220:244)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L260-L276" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks.HumannessTaskResult:to_sheets` (260:276)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/web/tasks.py#L287-L312" target="_blank" rel="noopener noreferrer">`biophi.humanization.web.tasks:humanness_task` (287:312)</a>


### Humanness Data Models
This component defines the core data structures used to represent and encapsulate humanness-related information for peptides, individual antibody chains, and complete antibodies. These models also include methods for data manipulation and formatting.


**Related Classes/Methods**:

- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L24-L34" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.OASisParams` (24:34)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L38-L56" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.PeptideHumanness` (38:56)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L60-L195" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness` (60:195)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L199-L255" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness` (199:255)</a>


### Core Humanness Analysis
Dedicated to calculating and assessing the humanness of antibody sequences, primarily using the OASis method. It includes statistical functions to determine humanness scores, identify non-human peptides, and generate related metrics, providing quantitative insights into sequence humanness.


**Related Classes/Methods**:

- `biophi.humanization.methods.humanness` (full file reference)
- `biophi.humanization.methods.stats` (full file reference)
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L273-L277" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.get_antibody_humanness` (273:277)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L319-L343" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.get_chain_humanness` (319:343)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L280-L316" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.get_chain_oasis_peptides` (280:316)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L346-L352" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.parse_peptide_humanness` (346:352)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L258-L270" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.chop_seq_peptides` (258:270)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L355-L367" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.get_oas_hits` (355:367)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L374-L377" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.get_fraction_subjects_bin` (374:377)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L76-L87" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_oasis_curve` (76:87)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L89-L90" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_oasis_identity` (89:90)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L92-L97" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_oasis_percentile` (92:97)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L107-L108" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_num_human_peptides` (107:108)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L110-L111" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_num_nonhuman_peptides` (110:111)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L113-L114" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.ChainHumanness.get_num_peptides` (113:114)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L203-L204" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness.get_oasis_identity` (203:204)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L206-L211" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness.get_oasis_percentile` (206:211)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L228-L230" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness.get_num_human_peptides` (228:230)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L232-L234" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness.get_num_nonhuman_peptides` (232:234)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanness.py#L225-L226" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanness.AntibodyHumanness.get_num_peptides` (225:226)</a>
- `biophi.humanization.methods.stats.get_oasis_percentile` (full file reference)
- `biophi.humanization.methods.stats.get_germline_family_residue_frequency` (full file reference)
- `biophi.humanization.methods.stats.get_chain_type_residue_frequency` (full file reference)
- `biophi.humanization.methods.stats._get_frequency` (full file reference)


### Humanization Core
This component is responsible for the actual humanization process of antibody sequences, distinct from humanness analysis. It contains algorithms for modifying sequences to increase their humanness.


**Related Classes/Methods**:

- `biophi.humanization.methods.humanization` (full file reference)
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanization.py#L42-L60" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanization.CDRGraftingHumanizationParams` (42:60)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanization.py#L64-L68" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanization.ManualHumanizationParams` (64:68)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanization.py#L24-L38" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanization.SapiensHumanizationParams` (24:38)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanization.py#L121-L155" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanization.AntibodyHumanization` (121:155)</a>
- <a href="https://github.com/Merck/BioPhi/blob/master/biophi/humanization/methods/humanization.py#L78-L117" target="_blank" rel="noopener noreferrer">`biophi.humanization.methods.humanization.ChainHumanization` (78:117)</a>




### [FAQ](https://github.com/CodeBoarding/GeneratedOnBoardings/tree/main?tab=readme-ov-file#faq)