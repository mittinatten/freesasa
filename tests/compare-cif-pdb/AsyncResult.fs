// from https://fsharpforfunandprofit.com/posts/elevated-world-5/
module AsyncResult

type AsyncResult<'a, 'e> = Async<Result<'a, 'e>>

let bind f xAsyncResult = async {
    let! xResult = xAsyncResult
    match xResult with
    | Ok x -> return! f x
    | Error err -> return (Error err)
    }
