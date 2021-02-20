module Download
open System
open System.IO
open System.Net


let downloadFile (url: Uri) path =
    async {
        if File.Exists(path) then
            return Ok ()
        else
            try
                let req = WebRequest.Create(url)
                use resp = req.GetResponse()
                use stream = resp.GetResponseStream()
                use fileStream = File.Create(path)
                do! stream.CopyToAsync(fileStream) |> Async.AwaitTask
                return Ok ()
            with ex -> return Error ("Error downloading " + url.ToString() + " " + ex.Message)
    }
