.pdxMI_data <- function()
{
  data(PDXMI)
  return(PDXMI)
}

createXevaModelClass <- function()
{
  PDXMI <- .pdxMI_data()
  pdxmiVar <- as.list(rep("character", length(PDXMI$id)))
  names(pdxmiVar) <- as.character(PDXMI$id)
  pdxmiPrototype <- as.list(rep(NA_character_, length(pdxmiVar)))
  names(pdxmiPrototype) <- names(pdxmiVar)
  ##----------------------------------------------------------------------------
  xevaVar <- list(model.id = "character", drug = "list", data="data.frame",
                  treatment.type = "character",
                  treatment.target="character")

  xevaPrototype <- list(treatment.type = NA_character_,
                        treatment.target=NA_character_)

  pdxmodSlote <- append(xevaVar, pdxmiVar)
  pdxmodPrototype <- append(xevaPrototype, pdxmiPrototype)

  PDXmodClass<-setClass("PDX_model",
                      slots=pdxmodSlote,
                      prototype=pdxmodPrototype
                      )

  return(PDXmodClass)
}


PDXmodClass <- createXevaModelClass()


## @export
#' @import methods
setMethod(f="show",
          signature="PDX_model",
          definition= function(object)
          {
            msg <- sprintf("model.id = %s\ndrug = %s\n",
                           slot(object, "model.id"),
                           slot(object, "drug")$join.name)
            cat(msg)
            cat(sprintf("\ndata\n"))
            print(head(slot(object, "data")))
            cat(sprintf("\n"))

            ##---------------
            otherSlot <- setdiff(slotNames(object), c("model.id", "drug", "data"))
            for(s in otherSlot)
            {
              if(!is.na(slot(object, s)))
              {
                cat(sprintf("%s = %s\n", s, slot(object, s)))
              }
            }
          })
