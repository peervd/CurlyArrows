import openai

def communicate_prompt(prompt,openai_key):
    # Initialize OpenAI client with API key
    client = openai.OpenAI(api_key=openai_key)

    # Call the API
    response = client.chat.completions.create(
        model='gpt-4o-mini',
        temperature = 0,
        messages=[{"role": "user", "content": prompt}])

    return response.choices[0].message.content

