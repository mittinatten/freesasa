module Main

open Command
open Download
open System
open System.IO

open FSharpPlus

type PdbCode = PdbCode of string

type FileFormat =
    | PDB
    | CIF

type PdbFile =
    { Path: string
      Format: FileFormat
      Code: PdbCode }

type Output =
    { Output: String
      Warnings: String
      ExitCode: int }

type CalculationError =
    | CalculationFailed of Output
    | DownloadFailed of String

type FreesasaResult = Result<Output, CalculationError>

type ComparisonError =
    | CalculationError of CalculationError list
    | CalculationDifference

type Comparison = (PdbCode * ComparisonError) option
type ResultComparer = PdbCode -> FreesasaResult -> FreesasaResult -> Comparison

let archiveUrl = Uri("https://files.rcsb.org/download/")
let dataDirectory = "./data/"
let freesasaProgram = "../../src/freesasa"

let getExtension format =
    match format with
    | PDB -> ".pdb"
    | CIF -> ".cif"

let downloadPdb code format =
    let getUrl (PdbCode code) =
        Uri(archiveUrl, code + (getExtension format))

    let getPath (PdbCode code) =
        dataDirectory + code + (getExtension format)

    async {
        let path = getPath code
        let url = getUrl code
        let! result = downloadFile url path

        return
            match result with
            | Ok () ->
                Ok
                    { Path = path
                      Format = format
                      Code = code }
            | Error e -> Error(DownloadFailed e)
    }

let normalizeCommandResult (commandResult: CommandResult) : Output =
    let stripSource output =
        String.split [ "\n" ] output
        |> Seq.filter (fun line -> not (line.StartsWith "source"))
        |> String.intercalate "\n"

    { Output = stripSource commandResult.StandardOutput
      Warnings = commandResult.StandardError.Replace(" ", "")
      ExitCode = commandResult.ExitCode }

let runPdbCalc args file =
    let log output =
        printf "\rCalculated %s (exit code %d)" file.Path output.ExitCode
        stdout.Flush()

    let allArgs =
        if file.Format = CIF then
            "--cif" :: args
        else
            args

    async {
        let! output =
            executeCommand freesasaProgram (file.Path :: allArgs)
            |> Async.map normalizeCommandResult

        do log output

        return
            match output with
            | { ExitCode = 0 } -> Ok output
            | _ -> Error(CalculationFailed output)
    }

let calcStructure code fileType args =
    downloadPdb code fileType
    |> AsyncResult.bind (runPdbCalc args)

let compareFormats args (compare: ResultComparer) code : Async<Comparison> =
    async {
        let! cifResult = calcStructure code CIF args
        let! pdbResult = calcStructure code PDB args
        return compare code cifResult pdbResult
    }

let compareResults : ResultComparer =
    fun code result1 result2 ->
        let err = CalculationDifference

        match (result1, result2) with
        | Ok v1, Ok v2 ->
            if (v1 = v2) then
                None
            else
                Some(code, err)
        | Error e1, Error e2 ->
            if (e1 = e2) then
                None
            else
                Some(code, err)
        | _ -> Some(code, err)

let checkArgs codes args =
    printf "Checking with args %A\n" args

    let computations =
        codes
        |> List.map (compareFormats args compareResults)

    Async.Parallel(computations, 8)
    |> Async.RunSynchronously
    |> Array.toList
    |> List.choose id
    |> fun errors ->
        if List.isEmpty errors then
            printfn "\r - No errors"
            None
        else
            printfn " - Comparison failed:\n %A\n" errors
            Some errors


let getCodes fileName =
    System.IO.File.ReadAllLines(fileName)
    |> Seq.map (String.trimEnd " " >> PdbCode)
    |> Seq.toList

[<EntryPoint>]
let main argv =
    if argv.Length < 1 then
        printfn "usage: pass input file as argument"
        1
    else
        if not (Directory.Exists dataDirectory) then
            do Directory.CreateDirectory(dataDirectory) |> ignore

        let codes = getCodes argv.[0]

        [ [] ] // ["--format=xml"] ["--separate-chains", "--separate-models"]
        |> List.map (checkArgs codes)
        |> List.choose id
        |> fun errList -> if List.isEmpty errList then 0 else 1
