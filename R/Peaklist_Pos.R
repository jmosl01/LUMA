#' @title Peaklist_Pos_db Tables
#'
#' @description Tables generated from the fish samples in positive ionization
#'   mode
#' @source
#'   https://https://github.com/jmosl01/LUMA/blob/develop/data-raw/Peaklist_Pos.SQLite
#'
#' @format list of tibbles \describe{
#'   \item{Annotated}{A tibble: 13 x 131, example for sum_features}
#'   \item{Combined_Isotopes_and_Adducts}{A tibble: 12 x 141, tests trim_minfrac
#'   method "monoMass"}
#'   \item{From_CAMERA}{A tibble: 23 x 120, tests trimming
#'   functions except trim_minfrac, tests calculation functions except
#'   calc_corrstat}
#'   \item{From_CAMERA_with_MinFrac}{A tibble: 23 x 121, tests
#'   trim_minfrac method "mz"}
#'   \item{Trimmed_by_CV}{A tibble: 3 x 142}
#'   \item{Trimmed_by_MinFrac}{A tibble: 2 x 142}
#'   \item{Trimmed_by_RT}{A tibble:
#'   13 x 121}
#'   \item{input_parsed}{A tibble: 13 x 138, tests calc_corrstat}
#'   \item{output_parsed}{A tibble: 12 x 139, tests sum_features}
#'    }
#' @examples
#' \dontrun{Peaklist_Pos}
"Peaklist_Pos"
