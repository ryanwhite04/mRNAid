import sendgrid
from os import getenv
from sendgrid.helpers.mail import Mail, Email, To, Content
# from dotenv import load_dotenv

# load_dotenv()
api_key = getenv('SENDGRID_API_KEY')
email_username = getenv('SENDGRID_EMAIL_USERNAME')
is_test_environment = getenv('FLASK_ENV') == 'testing'

sg = sendgrid.SendGridAPIClient(api_key)

def send_email(subject, body, recipient):
    if is_test_environment:
        print(f"[TEST ENVIRONMENT] Send email to {recipient} with subject '{subject}' and body '{body}'")
        return
    else:
        from_email = Email(email_username, "mRNA Server")  # the verified sender and name
        to_email = To(recipient) # recipients

        content = Content("text/plain", body)
        mail = Mail(from_email, to_email, subject, content).get()
        
        # Send an HTTP POST request to /mail/send
        response = sg.client.mail.send.post(request_body=mail)
        print(response.status_code)
        print(response.headers)
        print("Email sent successfully!")