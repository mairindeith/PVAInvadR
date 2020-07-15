# How to handle the decision inputs?

## Options:

### Which parameters to modify with controls?

- use `tk_select.list(choices, preselect = NULL, multiple = F, title = NULL)`
  - if selection = "", quit (no control params selected)

**From the Shiny app:**

```r
output$input.selectParams <- renderUI({
  checkboxGroupInput(inputId = "selectParams", label = "Select the parameters to vary between scenarios",
                     choices=c(paste0("Time-step to begin sampling pre-recruits, stanza ",1:input$nS,
                                      " (t.start.R.",1:input$nS,")"),
                               paste0("Time-step to begin sampling adults, gear ",1:input$n.gear,
                                      " (t.start.A.",1:input$n.gear,")"),
                               paste0("Time-step in the year that gear ",1:input$n.gear,
                                      " is fished (samp.A.",1:input$n.gear,")"),
                               paste0("Effort expended by each pre-recruit sampling gear, stanza ",1:input$nS,
                                      " (E.R.",1:input$nS,")"),
                               paste0("Proportion of pre-recruits removed with one unit of effort, stanza ",1:input$nS,
                                      " (U.R.",1:input$nS,")"),
                               paste0("Effort expended by each adult sampling gear, gear ",1:input$n.gear,
                                      " (E.A.",1:input$n.gear,")")})
```

**Frankensteined:**

```r
# input is the read-in parameters
# scen_names is just a character vector of scenario names

# Given input to decision_setup():
decision_setup = function(input, scen_names, list = T, csv = F, csv_path = NULL,
  gui = T, selected_params = NULL){
    # Identify all possible parameters that could be modified by scenarios

  # First, assess whether the GUI should pop up or not.
  all_params = c(
    paste0("t.start.R.", 1:input$nS),
    paste0("t.start.A.",1:input$n.gear),
    paste0("samp.A.",1:input$n.gear),
    paste0("E.R.",1:input$nS),
    paste0("U.R.",1:input$nS),
    paste0("E.A.",1:input$n.gear)
  )
  if(gui == F){
    # If the GUI does not pop up, selected parameters must be provided.
    if(is.null(selected_params)){
      stop("If gui = FALSE, which parameters to vary between scenarios must be provided as a vector of characters, selected_params.")
      return(NULL)
    } else { # if(!is.null(selected_params))
    # If selected_parameters are not NULL, are they in the set?
    # Check to see if the parameters in selected_params are valid:
    diffs = selected_params[!(which(selected_params %in% all_params))]
    if(!identical(diffs, character(0))){
      stop("Parameter names in selected_params not valid. decision_setup() failed.")
      return(NULL)
    }
  } # if()
  } else { # if(gui == T)
    param_desc=c(
      paste0("Time-step to begin sampling pre-recruits, stanza ",1:input$nS,
        " (t.start.R.",1:input$nS,")"),
      paste0("Time-step to begin sampling adults, gear ",1:input$n.gear,
        " (t.start.A.",1:input$n.gear,")"),
      paste0("Time-step in the year that gear ",1:input$n.gear,
        " is fished (samp.A.",1:input$n.gear,")"),
      paste0("Effort expended by each pre-recruit sampling gear, stanza ",
        1:input$nS," (E.R.",1:input$nS,")"),
      paste0("Proportion of pre-recruits removed with 1 unit of effort, stanza ",
        1:input$nS," (U.R.",1:input$nS,")"),
      paste0("Effort expended by each adult sampling gear, gear ",
        1:input$n.gear," (E.A.",1:input$n.gear,")")
    )
    selected_params = select.list(all_params, multiple = T,
      title = "Select the parameters to vary between scenarios")
    if(length(selected_params) == 0){
      stop("No parameters selected, decision_setup() failed.")
      return(NULL)
    }
  }
  if(list==T){
    scen_list = vector("list", length(scen_names)) # Outermost list, named after scenario names
    names(scen_list) = scen_names
    for(n in scen_names){
      scen_list[[n]] = vector("list", length(selected_params))
      names(scen_list[[n]]) = selected_params
    }
    return(scen_list)
  }
  if(csv==T){
    if(is.null(csv_file)){
      stop("No csv save path provided, decision_setup() failed.")
      return(NULL)
    }
    scen_df = data.frame(matrix(nrow=length(scen_names), ncol=length(selected_params)+1))
    colnames(scen_df) = c("Scenario", selected_params)
    scen_df$Scenario = scen_names
    write_csv(scen_df, csv_file)
  }
}

```

### CSV

- Save a template CSV to be filled in with emtpy columns:
  - R1C1 = blank
  - R1C2:Cn = scenario names (to be filled in by the user)
  - R2C1:Cn = parameter names
    - If the cells are left blank, the value should remain the default value
    - If the cell is filled, the value in that cell will be the new value used for decision-making

### Data frame

- `create_scenario()` function:

```r
create_scenarios = function(scen_names, parameters){
  scen_list = vector("list", length(scen_names)) # Outermost list, named after scenario names
  names(scen_list) = scen_names
  for(n in scen_names){
    scen_list[[n]] = vector("list", length(parameters))
    names(scen_list[[n]]) = parameters
  }
  return(scen_list)
}
```


  - Takes the user input parameter list as a set of names in a named list, the user can then fill in using, e.g., `scenA[["nT"]] = `
- Then unlist to treat in the same was as the CSV? Say, unlisting then putting into a dataframe which is then translated?
- After the selection pop-up, a dataframe is created with the same structure as the CSV, namely:
  - rownames = parameters (chosen with `tk_select.list()`)
  - columns can be defined in the argument as a character vector; this will create NAs

**Global argument to turn off GUI elements? How to carry over arguments when loading the library?**

**Setup a global "gui" variable in R that is binary, T/F, and read from there with each function?**


## May 21, 2020 - Notes

- Combine init and parameter input, maybe as an "input" object that is a named list?
  - Do this with 'init' function - put the parameters fed into init() as the first set in the resultant list
- Figure out:
  - decision():
    - how to auto-fill the NAs with trhe default values?
    - if not, then modify decision_setup()?

- How to convert from a nested named list into a dataframe:
  - Outer names (hierarchy 1) = scenario names, will become row names/column 1 data
  - Inner names (hierachy 2 under each 1) = parameter names and values
  - Could do a loop with rbind...

```{r}
# ...

scenNames = names(decision_list)
for(s in scenNames){
  n_param = length(decision_list[[s]])
  tmp_row = data.frame(matrix(nrow=1, ncol=n_param)) # May need to be n_param + 1 if the first column is the scenName
  lapply(rbind, decision_list ) #... ? )
  for(i in 1:length(decision_list[[s]])){
    param_name = name(decision_list[[s]][i])
    param_val = as.numeric(decision_list[[s]][i])
    tmp_row[,i] =
  }
}

# ...
```
  - c("scenName", )
