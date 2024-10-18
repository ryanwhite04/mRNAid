import sendgrid
from sendgrid.helpers.mail import Mail, Email, To, Content
from config import SENDGRID_API_KEY, SENDGRID_EMAIL_USERNAME

def send_email(subject, body, recipient):
    if SENDGRID_API_KEY == None:
        print("API_KEY not set, email notifications not enabled")
        return
    if SENDGRID_EMAIL_USERNAME == None:
        print("EMAIL_USERNAME not set, email notifications disabled")
        return
    sg = sendgrid.SendGridAPIClient(SENDGRID_API_KEY)
    from_email = Email(SENDGRID_EMAIL_USERNAME, "mRNA Server")  # the verified sender and name
    to_email = To(recipient)
    content = Content("text/plain", body)
    mail = Mail(from_email, to_email, subject, content).get()
    
    # Send an HTTP POST request to /mail/send
    return sg.client.mail.send.post(request_body=mail)