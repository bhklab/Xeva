#' @import utils
.pdxMI_data <- function()
{
  structure(list(Module = c("Clinical/patient", "Clinical/patient",
                            "Clinical/patient", "Clinical/patient", "Clinical/patient", "Clinical/patient",
                            "Clinical/patient", "Clinical/patient", "Clinical/patient", "Clinical/patient",
                            "Clinical/patient", "Clinical/tumor", "Clinical/tumor", "Clinical/tumor",
                            "Clinical/tumor", "Clinical/tumor", "Clinical/tumor", "Clinical/tumor",
                            "Clinical/tumor", "Clinical/tumor", "Clinical/tumor", "Clinical/tumor",
                            "Model creation", "Model creation", "Model creation", "Model creation",
                            "Model creation", "Model creation", "Model creation", "Model creation",
                            "Model creation", "Model quality assurance", "Model quality assurance",
                            "Model quality assurance", "Model quality assurance", "Model quality assurance",
                            "Model study", "Model study", "Model study", "Model study", "Model study",
                            "Model study", "Associated metadata", "Associated metadata",
                            "Associated metadata"), Field = c("Submitter patient ID", "Gender",
                                                              "Age", "Diagnosis", "Consent to share data", "Ethnicity/race",
                                                              "Current treatment drug", "Current treatment protocol (dose; details)",
                                                              "Prior treatment protocol", "Response to prior treatment", "Virology status",
                                                              "Submitter tumor ID", "Primary tumor tissue of origin", "Primary, metastasis, recurrence",
                                                              "Specimen tumor tissue", "Tissue histology", "Tumor grade; classification",
                                                              "Disease stage; classification", "Specific markers (diagnostic linked); platform",
                                                              "Is tumor from untreated patient?", "Original tumor sample type",
                                                              "Tumor from an existing PDX model? ID? Why sub-line?", "Submitter PDX ID",
                                                              "Mouse strain (and source)", "Strain immune system humanized?",
                                                              "Type of humanization", "Tumor preparation", "Injection type and site",
                                                              "Mouse treatment for engraftment", "Engraftment rate", "Engraftment time",
                                                              "Tumor characterization technology", "Tumor confirmed not to be of mouse/EBV origin",
                                                              "Response to standard of care (pharmacologic positive control)",
                                                              "Animal health status", "Passage QA performed", "Treatment, passage",
                                                              "Treatment protocol (dose; details)", "Treatment response", "Tumor OMICS: sample id; sample site; purity (mouse vs. human); technology; passage",
                                                              "Development of metastases in strain (Y/N, site); passage", "Lag time/doubling time of tumor",
                                                              "PDX model availability?", "Governance restriction for distribution",
                                                              "ID for associated publication, image, archived data (URL, PMID, DOI)"
                            ), Recommendation = c("Essential", "Essential", "Essential",
                                                  "Essential", "Essential", "Desirable", "Desirable", "Desirable",
                                                  "Desirable", "Desirable", "Desirable", "Essential", "Essential",
                                                  "Essential", "Essential", "Essential", "Essential", "Essential",
                                                  "Essential", "Essential", "Desirable", "Desirable", "Essential",
                                                  "Essential", "Essential", "Essential", "Essential", "Essential",
                                                  "Desirable", "Desirable", "Desirable", "Essential", "Essential",
                                                  "Desirable", "Desirable", "Essential", "Desirable", "Desirable",
                                                  "Desirable", "Desirable", "Desirable", "Desirable", "Desirable",
                                                  "Desirable", "Desirable"), id = c("patient.id", "patient.sex",
                                                                                    "patient.age", "patient.diagnosis", "patient.consent.to.share.data",
                                                                                    "patient.ethnicity", "patient.current.treatment.drug", "patient.current.treatment.protocol",
                                                                                    "patient.prior.treatment.protocol", "patient.response.to.prior.treatment",
                                                                                    "patient.virology.status", "tumor.id", "tumor.tissue.of.origin",
                                                                                    "tumor.primary.metastasis.recurrence", "tumor.specimen.tissue",
                                                                                    "tumor.tissue.histology", "tumor.tumor.grade", "tumor.disease.stage",
                                                                                    "tumor.specific.markers", "tumor.fom.untreated.patient", "tumor.original.sample.type",
                                                                                    "tumor.from.existing.pdx.model", "model.submitter.pdx.id", "model.mouse.strain.source",
                                                                                    "model.strain.immune.system.humanized", "model.type.of.humanization",
                                                                                    "model.tumor.preparation", "model.injection.type.and.site", "model.mouse.treatment.for.engraftment",
                                                                                    "model.engraftment.rate", "model.engraftment.time", "model.tumor.characterization.technology",
                                                                                    "model.tumor.confirmed.not.to.be.of.mouse.origin", "model.response.to.standard.of.care",
                                                                                    "model.animal.health.status", "model.passage.qa.performed", "model.treatment.passage",
                                                                                    "model.treatment.protocol", "model.treatment.response", "model.tumor.omics",
                                                                                    "model.development.of.metastases.in.strain", "model.doubling.time.of.tumor",
                                                                                    "pdx.model.availability", "governance.restriction.for.distribution",
                                                                                    "id.publication.data")), .Names = c("Module", "Field", "Recommendation",
                                                                                                                        "id"), row.names = c(NA, -45L), class = "data.frame")



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
                  meta="list",
                  treatment.type = "character",
                  treatment.target="character")

  xevaPrototype <- list(treatment.type = NA_character_,
                        treatment.target=NA_character_)

  pdxmodSlote <- append(xevaVar, pdxmiVar)
  pdxmodPrototype <- append(xevaPrototype, pdxmiPrototype)

  PDXmodClass<-setClass("pdxModel",
                      slots=pdxmodSlote,
                      prototype=pdxmodPrototype
                      )

  return(PDXmodClass)
}


PDXmodClass <- createXevaModelClass()


#' @import methods
setMethod(f="show",
          signature="pdxModel",
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
