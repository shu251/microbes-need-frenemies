# Factoring 

order_factor <- function(df){
  df |> 
  mutate(SAMPLE_NAME_ORDER = factor(raw_label, levels = order_sample_names, labels = label_sample_names),
         SAMPLE_NAME_ORDER_SHORT = factor(raw_label, levels = order_sample_names, labels = label_sample_names_SHORT),
         VENT_BIN_ORDER = factor(VENT_BIN, levels = order_bin),
         VENT_FIELD_ORDER = factor(VENT_FIELD, levels = order_field))
}
