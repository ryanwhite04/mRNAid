import sendgrid
from os import getenv
from sendgrid.helpers.mail import Mail, Email, To, Content
from dotenv import load_dotenv

load_dotenv()
api_key = getenv('SENDGRID_API_KEY')
email_username = getenv('SENDGRID_EMAIL_USERNAME')

sg = sendgrid.SendGridAPIClient(api_key)

def send_email(subject, body, recipient):
    from_email = Email(email_username, "mRNA Server")  # the verified sender and name
    to_email = To(recipient)  # recipients
    content = Content("text/plain", body)
    mail = Mail(from_email, to_email, subject, content).get()
    
    # Send an HTTP POST request to /mail/send
    response = sg.client.mail.send.post(request_body=mail)
    print(response.status_code)
    print(response.headers)
    print("Email sent successfully!")